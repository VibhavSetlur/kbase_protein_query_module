#!/usr/bin/env python
# -*- coding: utf-8 -*-
import datetime
import json
import os
import random as _random
import sys
import traceback
from getopt import getopt, GetoptError
from multiprocessing import Process
from os import environ
from wsgiref.simple_server import make_server

import requests as _requests
from jsonrpcbase import JSONRPCService, InvalidParamsError, KeywordError, \
    JSONRPCError, InvalidRequestError
from jsonrpcbase import ServerError as JSONServerError

from biokbase import log
from kbase_protein_query_module.authclient import KBaseAuth as _KBaseAuth

try:
    from ConfigParser import ConfigParser
except ImportError:
    from configparser import ConfigParser

DEPLOY = 'KB_DEPLOYMENT_CONFIG'
SERVICE = 'KB_SERVICE_NAME'
AUTH = 'auth-service-url'

# Note that the error fields do not match the 2.0 JSONRPC spec


def get_config_file():
    return environ.get(DEPLOY, None)


def get_service_name():
    return environ.get(SERVICE, None)


def get_config():
    if not get_config_file():
        return None
    retconfig = {}
    config = ConfigParser()
    config.read(get_config_file())
    for nameval in config.items(get_service_name() or 'kbase_protein_query_module'):
        retconfig[nameval[0]] = nameval[1]
    return retconfig

config = get_config()

from kbase_protein_query_module.kbase_protein_query_moduleImpl import kbase_protein_query_module  # noqa @IgnorePep8
impl_kbase_protein_query_module = kbase_protein_query_module(config)


class JSONObjectEncoder(json.JSONEncoder):

    def default(self, obj):
        if isinstance(obj, set):
            return list(obj)
        if isinstance(obj, frozenset):
            return list(obj)
        if hasattr(obj, 'toJSONable'):
            return obj.toJSONable()
        return json.JSONEncoder.default(self, obj)


class JSONRPCServiceCustom(JSONRPCService):

    def call(self, ctx, jsondata):
        """
        Calls jsonrpc service's method and returns its return value in a JSON
        string or None if there is none.

        Arguments:
        jsondata -- remote method call in jsonrpc format
        """
        result = self.call_py(ctx, jsondata)
        if result is not None:
            return json.dumps(result, cls=JSONObjectEncoder)

        return None

    def _call_method(self, ctx, request):
        """Calls given method with given params and returns it value."""
        method = self.method_data[request['method']]['method']
        params = request['params']
        result = None
        try:
            if isinstance(params, list):
                # Does it have enough arguments?
                if len(params) < self._man_args(method) - 1:
                    raise InvalidParamsError('not enough arguments')
                # Does it have too many arguments?
                if(not self._vargs(method) and len(params) >
                        self._max_args(method) - 1):
                    raise InvalidParamsError('too many arguments')

                result = method(ctx, *params)
            elif isinstance(params, dict):
                # Do not accept keyword arguments if the jsonrpc version is
                # not >=1.1.
                if request['jsonrpc'] < 11:
                    raise KeywordError

                result = method(ctx, **params)
            else:  # No params
                result = method(ctx)
        except JSONRPCError:
            raise
        except Exception as e:
            # log.exception('method %s threw an exception' % request['method'])
            # Exception was raised inside the method.
            newerr = JSONServerError()
            newerr.trace = traceback.format_exc()
            if len(e.args) == 1:
                newerr.data = repr(e.args[0])
            else:
                newerr.data = repr(e.args)
            raise newerr
        return result

    def call_py(self, ctx, jsondata):
        """
        Calls jsonrpc service's method and returns its return value in python
        object format or None if there is none.

        This method is same as call() except the return value is a python
        object instead of JSON string. This method is mainly only useful for
        debugging purposes.
        """
        rdata = jsondata
        # we already deserialize the json string earlier in the server code, no
        # need to do it again
#        try:
#            rdata = json.loads(jsondata)
#        except ValueError:
#            raise ParseError

        # set some default values for error handling
        request = self._get_default_vals()

        if isinstance(rdata, dict) and rdata:
            # It's a single request.
            self._fill_request(request, rdata)
            respond = self._handle_request(ctx, request)

            # Don't respond to notifications
            if respond is None:
                return None

            return respond
        elif isinstance(rdata, list) and rdata:
            # It's a batch.
            requests = []
            responds = []

            for rdata_ in rdata:
                # set some default values for error handling
                request_ = self._get_default_vals()
                self._fill_request(request_, rdata_)
                requests.append(request_)

            for request_ in requests:
                respond = self._handle_request(ctx, request_)
                # Don't respond to notifications
                if respond is not None:
                    responds.append(respond)

            if responds:
                return responds

            # Nothing to respond.
            return None
        else:
            # empty dict, list or wrong type
            raise InvalidRequestError

    def _handle_request(self, ctx, request):
        """Handles given request and returns its response."""
        if 'types' in self.method_data[request['method']]:
            self._validate_params_types(request['method'], request['params'])

        result = self._call_method(ctx, request)

        # Do not respond to notifications.
        if request['id'] is None:
            return None

        respond = {}
        self._fill_ver(request['jsonrpc'], respond)
        respond['result'] = result
        respond['id'] = request['id']

        return respond


class MethodContext(dict):

    def __init__(self, logger):
        self['client_ip'] = None
        self['user_id'] = None
        self['authenticated'] = None
        self['token'] = None
        self['module'] = None
        self['method'] = None
        self['call_id'] = None
        self['rpc_context'] = None
        self['provenance'] = None
        self._debug_levels = set([7, 8, 9, 'DEBUG', 'DEBUG2', 'DEBUG3'])
        self._logger = logger

    def log_err(self, message):
        self._log(log.ERR, message)

    def log_info(self, message):
        self._log(log.INFO, message)

    def log_debug(self, message, level=1):
        if level in self._debug_levels:
            pass
        else:
            level = int(level)
            if level < 1 or level > 3:
                raise ValueError("Illegal log level: " + str(level))
            level = level + 6
        self._log(level, message)

    def set_log_level(self, level):
        self._logger.set_log_level(level)

    def get_log_level(self):
        return self._logger.get_log_level()

    def clear_log_level(self):
        self._logger.clear_user_log_level()

    def _log(self, level, message):
        self._logger.log_message(level, message, self['client_ip'],
                                 self['user_id'], self['module'],
                                 self['method'], self['call_id'])

    def provenance(self):
        callbackURL = os.environ.get('SDK_CALLBACK_URL')
        if callbackURL:
            # OK, there's a callback server from which we can get provenance
            arg_hash = {'method': 'CallbackServer.get_provenance',
                        'params': [],
                        'version': '1.1',
                        'id': str(_random.random())[2:]
                        }
            body = json.dumps(arg_hash)
            response = _requests.post(callbackURL, data=body,
                                      timeout=60)
            response.encoding = 'utf-8'
            if response.status_code == 500:
                if ('content-type' in response.headers and
                        response.headers['content-type'] ==
                        'application/json'):
                    err = response.json()
                    if 'error' in err:
                        raise ServerError(**err['error'])
                    else:
                        raise ServerError('Unknown', 0, response.text)
                else:
                    raise ServerError('Unknown', 0, response.text)
            if not response.ok:
                response.raise_for_status()
            resp = response.json()
            if 'result' not in resp:
                raise ServerError('Unknown', 0,
                                  'An unknown server error occurred')
            return resp['result'][0]
        else:
            return self.get('provenance')


class ServerError(Exception):
    '''
    The call returned an error. Fields:
    name - the name of the error.
    code - the error code.
    message - a human readable error message.
    data - the server side stacktrace.
    '''

    def __init__(self, name, code, message, data=None, error=None):
        super(Exception, self).__init__(message)
        self.name = name
        self.code = code
        self.message = message if message else ''
        self.data = data or error or ''
        # data = JSON RPC 2.0, error = 1.1

    def __str__(self):
        return self.name + ': ' + str(self.code) + '. ' + self.message + \
            '\n' + self.data


def getIPAddress(environ):
    xFF = environ.get('HTTP_X_FORWARDED_FOR')
    realIP = environ.get('HTTP_X_REAL_IP')
    trustXHeaders = config is None or \
        config.get('dont_trust_x_ip_headers') != 'true'

    if (trustXHeaders):
        if (xFF):
            return xFF.split(',')[0].strip()
        if (realIP):
            return realIP.strip()
    return environ.get('REMOTE_ADDR')


class Application(object):
    # Wrap the wsgi handler in a class definition so that we can
    # do some initialization and avoid regenerating stuff over
    # and over

    def logcallback(self):
        self.serverlog.set_log_file(self.userlog.get_log_file())

    def log(self, level, context, message):
        self.serverlog.log_message(level, message, context['client_ip'],
                                   context['user_id'], context['module'],
                                   context['method'], context['call_id'])

    def __init__(self):
        submod = get_service_name() or 'kbase_protein_query_module'
        self.userlog = log.log(
            submod, ip_address=True, authuser=True, module=True, method=True,
            call_id=True, changecallback=self.logcallback,
            config=get_config_file())
        self.serverlog = log.log(
            submod, ip_address=True, authuser=True, module=True, method=True,
            call_id=True, logfile=self.userlog.get_log_file())
        self.serverlog.set_log_level(6)
        self.rpc_service = JSONRPCServiceCustom()
        self.method_authentication = dict()
        self.rpc_service.add(impl_kbase_protein_query_module.run_kbase_protein_query_module,
                             name='kbase_protein_query_module.run_kbase_protein_query_module',
                             types=[dict])
        self.method_authentication['kbase_protein_query_module.run_kbase_protein_query_module'] = 'required'  # noqa
        self.rpc_service.add(impl_kbase_protein_query_module.check_protein_existence,
                             name='kbase_protein_query_module.check_protein_existence',
                             types=[dict])
        self.method_authentication['kbase_protein_query_module.check_protein_existence'] = 'required'  # noqa
        self.rpc_service.add(impl_kbase_protein_query_module.generate_protein_embedding,
                             name='kbase_protein_query_module.generate_protein_embedding',
                             types=[dict])
        self.method_authentication['kbase_protein_query_module.generate_protein_embedding'] = 'required'  # noqa
        self.rpc_service.add(impl_kbase_protein_query_module.assign_family_fast,
                             name='kbase_protein_query_module.assign_family_fast',
                             types=[dict])
        self.method_authentication['kbase_protein_query_module.assign_family_fast'] = 'required'  # noqa
        self.rpc_service.add(impl_kbase_protein_query_module.find_top_matches_from_embedding,
                             name='kbase_protein_query_module.find_top_matches_from_embedding',
                             types=[dict])
        self.method_authentication['kbase_protein_query_module.find_top_matches_from_embedding'] = 'required'  # noqa
        self.rpc_service.add(impl_kbase_protein_query_module.summarize_and_visualize_results,
                             name='kbase_protein_query_module.summarize_and_visualize_results',
                             types=[dict])
        self.method_authentication['kbase_protein_query_module.summarize_and_visualize_results'] = 'required'  # noqa
        self.rpc_service.add(impl_kbase_protein_query_module.status,
                             name='kbase_protein_query_module.status',
                             types=[dict])
        authurl = config.get(AUTH) if config else None
        self.auth_client = _KBaseAuth(authurl)

    def __call__(self, environ, start_response):
        # Context object, equivalent to the perl impl CallContext
        ctx = MethodContext(self.userlog)
        ctx['client_ip'] = getIPAddress(environ)
        status = '500 Internal Server Error'

        try:
            body_size = int(environ.get('CONTENT_LENGTH', 0))
        except (ValueError):
            body_size = 0
        if environ['REQUEST_METHOD'] == 'OPTIONS':
            # we basically do nothing and just return headers
            status = '200 OK'
            rpc_result = ""
        else:
            request_body = environ['wsgi.input'].read(body_size)
            try:
                req = json.loads(request_body)
            except ValueError as ve:
                err = {'error': {'code': -32700,
                                 'name': "Parse error",
                                 'message': str(ve),
                                 }
                       }
                rpc_result = self.process_error(err, ctx, {'version': '1.1'})
            else:
                ctx['module'], ctx['method'] = req['method'].split('.')
                ctx['call_id'] = req['id']
                ctx['rpc_context'] = {
                    'call_stack': [{'time': self.now_in_utc(),
                                    'method': req['method']}
                                   ]
                }
                prov_action = {'service': ctx['module'],
                               'method': ctx['method'],
                               'method_params': req['params']
                               }
                ctx['provenance'] = [prov_action]
                try:
                    token = environ.get('HTTP_AUTHORIZATION')
                    # parse out the method being requested and check if it
                    # has an authentication requirement
                    method_name = req['method']
                    auth_req = self.method_authentication.get(
                        method_name, 'none')
                    if auth_req != 'none':
                        if token is None and auth_req == 'required':
                            err = JSONServerError()
                            err.data = (
                                'Authentication required for ' +
                                'kbase_protein_query_module ' +
                                'but no authentication header was passed')
                            raise err
                        elif token is None and auth_req == 'optional':
                            pass
                        else:
                            try:
                                user = self.auth_client.get_user(token)
                                ctx['user_id'] = user
                                ctx['authenticated'] = 1
                                ctx['token'] = token
                            except Exception as e:
                                if auth_req == 'required':
                                    err = JSONServerError()
                                    err.data = \
                                        "Token validation failed: %s" % e
                                    raise err
                    if (environ.get('HTTP_X_FORWARDED_FOR')):
                        self.log(log.INFO, ctx, 'X-Forwarded-For: ' +
                                 environ.get('HTTP_X_FORWARDED_FOR'))
                    self.log(log.INFO, ctx, 'start method')
                    rpc_result = self.rpc_service.call(ctx, req)
                    self.log(log.INFO, ctx, 'end method')
                    status = '200 OK'
                except JSONRPCError as jre:
                    err = {'error': {'code': jre.code,
                                     'name': jre.message,
                                     'message': jre.data
                                     }
                           }
                    trace = jre.trace if hasattr(jre, 'trace') else None
                    rpc_result = self.process_error(err, ctx, req, trace)
                except Exception:
                    err = {'error': {'code': 0,
                                     'name': 'Unexpected Server Error',
                                     'message': 'An unexpected server error ' +
                                                'occurred',
                                     }
                           }
                    rpc_result = self.process_error(err, ctx, req,
                                                    traceback.format_exc())

        # print('Request method was %s\n' % environ['REQUEST_METHOD'])
        # print('Environment dictionary is:\n%s\n' % pprint.pformat(environ))
        # print('Request body was: %s' % request_body)
        # print('Result from the method call is:\n%s\n' % \
        #    pprint.pformat(rpc_result))

        if rpc_result:
            response_body = rpc_result
        else:
            response_body = ''

        response_headers = [
            ('Access-Control-Allow-Origin', '*'),
            ('Access-Control-Allow-Headers', environ.get(
                'HTTP_ACCESS_CONTROL_REQUEST_HEADERS', 'authorization')),
            ('content-type', 'application/json'),
            ('content-length', str(len(response_body)))]
        start_response(status, response_headers)
        return [response_body.encode('utf8')]

    def process_error(self, error, context, request, trace=None):
        if trace:
            self.log(log.ERR, context, trace.split('\n')[0:-1])
        if 'id' in request:
            error['id'] = request['id']
        if 'version' in request:
            error['version'] = request['version']
            e = error['error'].get('error')
            if not e:
                error['error']['error'] = trace
        elif 'jsonrpc' in request:
            error['jsonrpc'] = request['jsonrpc']
            error['error']['data'] = trace
        else:
            error['version'] = '1.0'
            error['error']['error'] = trace
        return json.dumps(error)

    def now_in_utc(self):
        # noqa Taken from http://stackoverflow.com/questions/3401428/how-to-get-an-isoformat-datetime-string-including-the-default-timezone @IgnorePep8
        dtnow = datetime.datetime.now()
        dtutcnow = datetime.datetime.utcnow()
        delta = dtnow - dtutcnow
        hh, mm = divmod((delta.days * 24 * 60 * 60 + delta.seconds + 30) // 60,
                        60)
        return "%s%+02d:%02d" % (dtnow.isoformat(), hh, mm)

application = Application()

# This is the uwsgi application dictionary. On startup uwsgi will look
# for this dict and pull its configuration from here.
# This simply lists where to "mount" the application in the URL path
#
# This uwsgi module "magically" appears when running the app within
# uwsgi and is not available otherwise, so wrap an exception handler
# around it
#
# To run this server in uwsgi with 4 workers listening on port 9999 use:
# uwsgi -M -p 4 --http :9999 --wsgi-file _this_file_
# To run a using the single threaded python BaseHTTP service
# listening on port 9999 by default execute this file
#
try:
    import uwsgi
# Before we do anything with the application, see if the
# configs specify patching all std routines to be asynch
# *ONLY* use this if you are going to wrap the service in
# a wsgi container that has enabled gevent, such as
# uwsgi with the --gevent option
    if config is not None and config.get('gevent_monkeypatch_all', False):
        print("Monkeypatching std libraries for async")
        from gevent import monkey
        monkey.patch_all()
    uwsgi.applications = {'': application}
except ImportError:
    # Not available outside of wsgi, ignore
    pass

_proc = None


def start_server(host='localhost', port=0, newprocess=False):
    '''
    By default, will start the server on localhost on a system assigned port
    in the main thread. Excecution of the main thread will stay in the server
    main loop until interrupted. To run the server in a separate process, and
    thus allow the stop_server method to be called, set newprocess = True. This
    will also allow returning of the port number.'''

    global _proc
    if _proc:
        raise RuntimeError('server is already running')
    httpd = make_server(host, port, application)
    port = httpd.server_address[1]
    print("Listening on port %s" % port)
    if newprocess:
        _proc = Process(target=httpd.serve_forever)
        _proc.daemon = True
        _proc.start()
    else:
        httpd.serve_forever()
    return port


def stop_server():
    global _proc
    _proc.terminate()
    _proc = None


def process_async_cli(input_file_path, output_file_path, token):
    exit_code = 0
    with open(input_file_path) as data_file:
        req = json.load(data_file)
    if 'version' not in req:
        req['version'] = '1.1'
    if 'id' not in req:
        req['id'] = str(_random.random())[2:]
    ctx = MethodContext(application.userlog)
    if token:
        user = application.auth_client.get_user(token)
        ctx['user_id'] = user
        ctx['authenticated'] = 1
        ctx['token'] = token
    if 'context' in req:
        ctx['rpc_context'] = req['context']
    ctx['CLI'] = 1
    ctx['module'], ctx['method'] = req['method'].split('.')
    prov_action = {'service': ctx['module'], 'method': ctx['method'],
                   'method_params': req['params']}
    ctx['provenance'] = [prov_action]
    resp = None
    try:
        resp = application.rpc_service.call_py(ctx, req)
    except JSONRPCError as jre:
        trace = jre.trace if hasattr(jre, 'trace') else None
        resp = {'id': req['id'],
                'version': req['version'],
                'error': {'code': jre.code,
                          'name': jre.message,
                          'message': jre.data,
                          'error': trace}
                }
    except Exception:
        trace = traceback.format_exc()
        resp = {'id': req['id'],
                'version': req['version'],
                'error': {'code': 0,
                          'name': 'Unexpected Server Error',
                          'message': 'An unexpected server error occurred',
                          'error': trace}
                }
    if 'error' in resp:
        exit_code = 500
    with open(output_file_path, "w") as f:
        f.write(json.dumps(resp, cls=JSONObjectEncoder))
    return exit_code

if __name__ == "__main__":
    if (len(sys.argv) >= 3 and len(sys.argv) <= 4 and
            os.path.isfile(sys.argv[1])):
        token = None
        if len(sys.argv) == 4:
            if os.path.isfile(sys.argv[3]):
                with open(sys.argv[3]) as token_file:
                    token = token_file.read()
            else:
                token = sys.argv[3]
        sys.exit(process_async_cli(sys.argv[1], sys.argv[2], token))
    try:
        opts, args = getopt(sys.argv[1:], "", ["port=", "host="])
    except GetoptError as err:
        # print help information and exit:
        print(str(err))  # will print something like "option -a not recognized"
        sys.exit(2)
    port = 9999
    host = 'localhost'
    for o, a in opts:
        if o == '--port':
            port = int(a)
        elif o == '--host':
            host = a
            print("Host set to %s" % host)
        else:
            assert False, "unhandled option"

    start_server(host=host, port=port)
#    print("Listening on port %s" % port)
#    httpd = make_server( host, port, application)
#
#    httpd.serve_forever()
