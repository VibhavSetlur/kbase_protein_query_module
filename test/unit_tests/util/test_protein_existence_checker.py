"""
Unit tests for ProteinExistenceChecker
"""
import pytest
from unittest.mock import Mock, patch
from lib.kbase_protein_query_module.src.util.storage.protein_existence_checker import ProteinExistenceChecker


class TestProteinExistenceChecker:
    """Test cases for ProteinExistenceChecker"""
    
    def test_existence_checker_initialization(self):
        """Test ProteinExistenceChecker initializes correctly"""
        checker = ProteinExistenceChecker()
        assert checker.cache == {}
        assert checker.cache_ttl == 3600  # Default TTL
    
    def test_existence_checker_custom_ttl(self):
        """Test ProteinExistenceChecker with custom TTL"""
        checker = ProteinExistenceChecker(cache_ttl=1800)
        assert checker.cache_ttl == 1800
    
    @patch('requests.get')
    def test_check_protein_exists_uniprot(self, mock_get):
        """Test checking protein existence via UniProt"""
        checker = ProteinExistenceChecker()
        
        # Mock successful response
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {
            "results": [{"primaryAccession": "P12345"}]
        }
        mock_get.return_value = mock_response
        
        exists, metadata = checker.check_protein_exists("P12345")
        
        assert exists is True
        assert metadata["source"] == "uniprot"
        assert metadata["accession"] == "P12345"
    
    @patch('requests.get')
    def test_check_protein_not_exists_uniprot(self, mock_get):
        """Test checking non-existent protein via UniProt"""
        checker = ProteinExistenceChecker()
        
        # Mock not found response
        mock_response = Mock()
        mock_response.status_code = 404
        mock_get.return_value = mock_response
        
        exists, metadata = checker.check_protein_exists("P99999")
        
        assert exists is False
        assert metadata["source"] == "uniprot"
        assert metadata["error"] == "not_found"
    
    @patch('requests.get')
    def test_check_protein_exists_api_error(self, mock_get):
        """Test handling API errors"""
        checker = ProteinExistenceChecker()
        
        # Mock API error
        mock_get.side_effect = Exception("API Error")
        
        exists, metadata = checker.check_protein_exists("P12345")
        
        assert exists is False
        assert metadata["source"] == "uniprot"
        assert "error" in metadata
    
    def test_check_protein_exists_cached(self):
        """Test checking cached protein existence"""
        checker = ProteinExistenceChecker()
        
        # Add to cache
        checker.cache["P12345"] = {
            "exists": True,
            "metadata": {"source": "uniprot", "accession": "P12345"},
            "timestamp": checker._get_timestamp()
        }
        
        exists, metadata = checker.check_protein_exists("P12345")
        
        assert exists is True
        assert metadata["source"] == "uniprot"
        assert metadata["accession"] == "P12345"
    
    def test_check_protein_exists_cache_expired(self):
        """Test checking protein with expired cache"""
        checker = ProteinExistenceChecker(cache_ttl=1)  # 1 second TTL
        
        # Add expired entry to cache
        checker.cache["P12345"] = {
            "exists": True,
            "metadata": {"source": "uniprot", "accession": "P12345"},
            "timestamp": checker._get_timestamp() - 2  # 2 seconds ago
        }
        
        with patch('requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.json.return_value = {
                "results": [{"primaryAccession": "P12345"}]
            }
            mock_get.return_value = mock_response
            
            exists, metadata = checker.check_protein_exists("P12345")
            
            assert exists is True
            mock_get.assert_called_once()  # Should make API call
    
    def test_check_multiple_proteins(self):
        """Test checking multiple proteins"""
        checker = ProteinExistenceChecker()
        
        with patch('requests.get') as mock_get:
            mock_response = Mock()
            mock_response.status_code = 200
            mock_response.json.return_value = {"entryType": "UniProtKB reviewed (Swiss-Prot)"}
            mock_get.return_value = mock_response
            
            results = checker.check_multiple_proteins(["P12345", "P67890"])
            
            assert len(results) == 2
            assert "P12345" in results
            assert "P67890" in results
            assert results["P12345"]["exists"] is True
            assert results["P67890"]["exists"] is True
    
    def test_clear_cache(self):
        """Test clearing the cache"""
        checker = ProteinExistenceChecker()
        
        # Add some entries to cache
        checker.cache["P12345"] = {"exists": True, "metadata": {}}
        checker.cache["P67890"] = {"exists": False, "metadata": {}}
        
        assert len(checker.cache) == 2
        
        checker.clear_cache()
        
        assert len(checker.cache) == 0
    
    def test_get_cache_stats(self):
        """Test getting cache statistics"""
        checker = ProteinExistenceChecker()
        
        # Add some entries to cache
        checker.cache["P12345"] = {"exists": True, "metadata": {}}
        checker.cache["P67890"] = {"exists": False, "metadata": {}}
        
        stats = checker.get_cache_stats()
        
        assert stats["total_entries"] == 2
        assert stats["cache_hit_rate"] == 0.0  # No hits yet
        assert stats["cache_size_bytes"] > 0
    
    def test_validate_protein_id(self):
        """Test validating protein ID format"""
        checker = ProteinExistenceChecker()
        
        # Valid UniProt IDs
        assert checker._validate_protein_id("P12345") is True
        assert checker._validate_protein_id("Q67890") is True
        assert checker._validate_protein_id("A1B2C3") is True
        
        # Invalid IDs
        assert checker._validate_protein_id("invalid") is False
        assert checker._validate_protein_id("12345") is False
        assert checker._validate_protein_id("") is False
    
    def test_get_timestamp(self):
        """Test getting current timestamp"""
        checker = ProteinExistenceChecker()
        
        timestamp = checker._get_timestamp()
        assert isinstance(timestamp, (int, float))
        assert timestamp > 0
    
    def test_is_cache_entry_valid(self):
        """Test checking if cache entry is valid"""
        checker = ProteinExistenceChecker(cache_ttl=3600)
        
        current_time = checker._get_timestamp()
        
        # Valid entry
        entry = {"timestamp": current_time}
        assert checker._is_cache_entry_valid(entry) is True
        
        # Expired entry
        entry = {"timestamp": current_time - 7200}  # 2 hours ago
        assert checker._is_cache_entry_valid(entry) is False
    
    def test_batch_check_with_mixed_results(self):
        """Test batch checking with mixed results"""
        checker = ProteinExistenceChecker()
        
        with patch('requests.get') as mock_get:
            # Mock different responses for different proteins
            def mock_get_side_effect(url, **kwargs):
                mock_response = Mock()
                if "P12345" in url:
                    mock_response.status_code = 200
                    mock_response.json.return_value = {"entryType": "UniProtKB reviewed (Swiss-Prot)"}
                else:
                    mock_response.status_code = 404
                return mock_response
            
            mock_get.side_effect = mock_get_side_effect
            
            results = checker.check_multiple_proteins(["P12345", "P99999"])
            
            assert results["P12345"]["exists"] is True
            assert results["P99999"]["exists"] is False
