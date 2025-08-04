package us.kbase.kbaseproteinquerymodule;

import java.util.List;
import java.util.Map;

/**
 * <p>Original spec-file type: ReportResults</p>
 * <pre>
 * ReportResults is a reference to a hash where the following keys are defined:
 * report_name has a value which is a string
 * report_ref has a value which is a string
 * </pre>
 */
public class ReportResults {
    private String reportName;
    private String reportRef;

    public ReportResults() {
    }

    public ReportResults(String reportName, String reportRef) {
        this.reportName = reportName;
        this.reportRef = reportRef;
    }

    public String getReportName() {
        return reportName;
    }

    public void setReportName(String reportName) {
        this.reportName = reportName;
    }

    public String getReportRef() {
        return reportRef;
    }

    public void setReportRef(String reportRef) {
        this.reportRef = reportRef;
    }
} 