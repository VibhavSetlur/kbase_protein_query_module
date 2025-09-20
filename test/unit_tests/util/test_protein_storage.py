"""
Unit tests for ProteinStorage
"""
import pytest
import numpy as np
from unittest.mock import Mock, patch, mock_open
from lib.kbase_protein_query_module.src.util.storage.protein_storage import ProteinStorage


class TestProteinStorage:
    """Test cases for ProteinStorage"""
    
    def test_protein_storage_initialization(self):
        """Test ProteinStorage initializes correctly"""
        storage = ProteinStorage()
        assert storage.proteins == {}
        assert storage.embeddings == {}
        assert storage.metadata == {}
        assert storage.index is None
    
    def test_add_protein(self):
        """Test adding a protein"""
        storage = ProteinStorage()
        
        protein_id = "P12345"
        sequence = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        metadata = {"source": "uniprot", "length": len(sequence)}
        
        storage.add_protein(protein_id, sequence, metadata)
        
        assert protein_id in storage.proteins
        assert storage.proteins[protein_id] == sequence
        assert storage.metadata[protein_id] == metadata
    
    def test_add_protein_with_embedding(self):
        """Test adding a protein with embedding"""
        storage = ProteinStorage()
        
        protein_id = "P12345"
        sequence = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        embedding = np.random.rand(128)
        metadata = {"source": "uniprot"}
        
        storage.add_protein(protein_id, sequence, metadata, embedding)
        
        assert protein_id in storage.proteins
        assert protein_id in storage.embeddings
        assert np.array_equal(storage.embeddings[protein_id], embedding)
    
    def test_get_protein(self):
        """Test getting a protein"""
        storage = ProteinStorage()
        
        protein_id = "P12345"
        sequence = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        
        storage.add_protein(protein_id, sequence)
        
        retrieved = storage.get_protein(protein_id)
        assert retrieved == sequence
        
        # Test non-existent protein
        retrieved = storage.get_protein("non_existent")
        assert retrieved is None
    
    def test_get_embedding(self):
        """Test getting a protein embedding"""
        storage = ProteinStorage()
        
        protein_id = "P12345"
        embedding = np.random.rand(128)
        
        storage.embeddings[protein_id] = embedding
        
        retrieved = storage.get_embedding(protein_id)
        assert np.array_equal(retrieved, embedding)
        
        # Test non-existent embedding
        retrieved = storage.get_embedding("non_existent")
        assert retrieved is None
    
    def test_get_metadata(self):
        """Test getting protein metadata"""
        storage = ProteinStorage()
        
        protein_id = "P12345"
        metadata = {"source": "uniprot", "length": 100}
        
        storage.metadata[protein_id] = metadata
        
        retrieved = storage.get_metadata(protein_id)
        assert retrieved == metadata
        
        # Test non-existent metadata
        retrieved = storage.get_metadata("non_existent")
        assert retrieved is None
    
    def test_list_proteins(self):
        """Test listing all proteins"""
        storage = ProteinStorage()
        
        storage.add_protein("P12345", "sequence1")
        storage.add_protein("P67890", "sequence2")
        
        proteins = storage.list_proteins()
        assert len(proteins) == 2
        assert "P12345" in proteins
        assert "P67890" in proteins
    
    def test_remove_protein(self):
        """Test removing a protein"""
        storage = ProteinStorage()
        
        protein_id = "P12345"
        sequence = "MKTVRQERLKSIVRILERSKEPVSGAQLAEELSVSRQVIVQDIAYLRSLGYNIVATPRGYVLAGG"
        embedding = np.random.rand(128)
        metadata = {"source": "uniprot"}
        
        storage.add_protein(protein_id, sequence, metadata, embedding)
        
        assert protein_id in storage.proteins
        assert protein_id in storage.embeddings
        assert protein_id in storage.metadata
        
        storage.remove_protein(protein_id)
        
        assert protein_id not in storage.proteins
        assert protein_id not in storage.embeddings
        assert protein_id not in storage.metadata
    
    def test_clear_all(self):
        """Test clearing all data"""
        storage = ProteinStorage()
        
        storage.add_protein("P12345", "sequence1")
        storage.add_protein("P67890", "sequence2")
        
        assert len(storage.proteins) == 2
        
        storage.clear_all()
        
        assert len(storage.proteins) == 0
        assert len(storage.embeddings) == 0
        assert len(storage.metadata) == 0
    
    def test_get_stats(self):
        """Test getting storage statistics"""
        storage = ProteinStorage()
        
        storage.add_protein("P12345", "sequence1", {"source": "uniprot"})
        storage.add_protein("P67890", "sequence2", {"source": "uniprot"}, np.random.rand(128))
        
        stats = storage.get_stats()
        
        assert stats["total_proteins"] == 2
        assert stats["proteins_with_embeddings"] == 1
        assert stats["proteins_with_metadata"] == 2
    
    def test_build_index(self):
        """Test building similarity index"""
        storage = ProteinStorage()
        
        # Add some proteins with embeddings
        for i in range(5):
            protein_id = f"P{i:05d}"
            embedding = np.random.rand(128)
            storage.add_protein(protein_id, f"sequence{i}", embedding=embedding)
        
        storage.build_index()
        
        assert storage.index is not None
    
    def test_search_similar(self):
        """Test searching for similar proteins"""
        storage = ProteinStorage()
        
        # Add some proteins with embeddings
        for i in range(5):
            protein_id = f"P{i:05d}"
            embedding = np.random.rand(128)
            storage.add_protein(protein_id, f"sequence{i}", embedding=embedding)
        
        storage.build_index()
        
        query_embedding = np.random.rand(128)
        results = storage.search_similar(query_embedding, k=3)
        
        assert len(results) <= 3
        assert all(isinstance(result, tuple) for result in results)
        assert all(len(result) == 2 for result in results)  # (protein_id, similarity)
    
    def test_save_to_file(self):
        """Test saving storage to file"""
        storage = ProteinStorage()
        
        storage.add_protein("P12345", "sequence1", {"source": "uniprot"})
        
        with patch('builtins.open', mock_open()) as mock_file:
            with patch('pickle.dump', Mock()) as mock_dump:
                storage.save_to_file("test_storage.pkl")
                mock_dump.assert_called_once()
    
    def test_load_from_file(self):
        """Test loading storage from file"""
        storage = ProteinStorage()
        
        mock_data = {
            "proteins": {"P12345": "sequence1"},
            "embeddings": {},
            "metadata": {"P12345": {"source": "uniprot"}}
        }
        
        with patch('builtins.open', mock_open()) as mock_file:
            with patch('pickle.load', return_value=mock_data) as mock_load:
                storage.load_from_file("test_storage.pkl")
                assert "P12345" in storage.proteins
                assert storage.proteins["P12345"] == "sequence1"
