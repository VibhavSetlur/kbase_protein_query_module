"""
Unit tests for ProteinFamilyAssigner
"""
import pytest
import numpy as np
from unittest.mock import Mock, patch
from lib.kbase_protein_query_module.src.util.storage.protein_family_assigner import ProteinFamilyAssigner


class TestProteinFamilyAssigner:
    """Test cases for ProteinFamilyAssigner"""
    
    def test_family_assigner_initialization(self):
        """Test ProteinFamilyAssigner initializes correctly"""
        assigner = ProteinFamilyAssigner()
        assert assigner.families == {}
        assert assigner.centroids == {}
        assert assigner.threshold == 0.7  # Default threshold
    
    def test_family_assigner_custom_threshold(self):
        """Test ProteinFamilyAssigner with custom threshold"""
        assigner = ProteinFamilyAssigner(threshold=0.8)
        assert assigner.threshold == 0.8
    
    def test_create_family(self):
        """Test creating a new family"""
        assigner = ProteinFamilyAssigner()
        
        family_id = "family_001"
        centroid = np.random.rand(128)
        
        assigner.create_family(family_id, centroid)
        
        assert family_id in assigner.families
        assert family_id in assigner.centroids
        assert np.array_equal(assigner.centroids[family_id], centroid)
    
    def test_add_protein_to_family(self):
        """Test adding a protein to a family"""
        assigner = ProteinFamilyAssigner()
        
        family_id = "family_001"
        protein_id = "P12345"
        embedding = np.random.rand(128)
        
        assigner.create_family(family_id, np.random.rand(128))
        assigner.add_protein_to_family(family_id, protein_id, embedding)
        
        assert protein_id in assigner.families[family_id]
        assert np.array_equal(assigner.families[family_id][protein_id], embedding)
    
    def test_assign_protein_to_family(self):
        """Test assigning a protein to the best matching family"""
        assigner = ProteinFamilyAssigner(threshold=0.5)
        
        # Create some families
        family1_centroid = np.array([1.0, 0.0, 0.0, 0.0])
        family2_centroid = np.array([0.0, 1.0, 0.0, 0.0])
        
        assigner.create_family("family_001", family1_centroid)
        assigner.create_family("family_002", family2_centroid)
        
        # Test protein that matches family1
        protein_embedding = np.array([0.9, 0.1, 0.0, 0.0])
        family_id = assigner.assign_protein_to_family("P12345", protein_embedding)
        
        assert family_id == "family_001"
        assert "P12345" in assigner.families["family_001"]
    
    def test_assign_protein_no_match(self):
        """Test assigning a protein when no family matches threshold"""
        assigner = ProteinFamilyAssigner(threshold=0.9)  # High threshold
        
        # Create a family
        family_centroid = np.array([1.0, 0.0, 0.0, 0.0])
        assigner.create_family("family_001", family_centroid)
        
        # Test protein that doesn't match well
        protein_embedding = np.array([0.0, 1.0, 0.0, 0.0])
        family_id = assigner.assign_protein_to_family("P12345", protein_embedding)
        
        assert family_id is None
        assert "P12345" not in assigner.families["family_001"]
    
    def test_get_family_proteins(self):
        """Test getting proteins in a family"""
        assigner = ProteinFamilyAssigner()
        
        family_id = "family_001"
        assigner.create_family(family_id, np.random.rand(128))
        
        assigner.add_protein_to_family(family_id, "P12345", np.random.rand(128))
        assigner.add_protein_to_family(family_id, "P67890", np.random.rand(128))
        
        proteins = assigner.get_family_proteins(family_id)
        assert len(proteins) == 2
        assert "P12345" in proteins
        assert "P67890" in proteins
    
    def test_get_family_proteins_nonexistent(self):
        """Test getting proteins from non-existent family"""
        assigner = ProteinFamilyAssigner()
        
        proteins = assigner.get_family_proteins("non_existent")
        assert proteins is None
    
    def test_get_family_centroid(self):
        """Test getting family centroid"""
        assigner = ProteinFamilyAssigner()
        
        family_id = "family_001"
        centroid = np.random.rand(128)
        assigner.create_family(family_id, centroid)
        
        retrieved_centroid = assigner.get_family_centroid(family_id)
        assert np.array_equal(retrieved_centroid, centroid)
    
    def test_get_family_centroid_nonexistent(self):
        """Test getting centroid from non-existent family"""
        assigner = ProteinFamilyAssigner()
        
        centroid = assigner.get_family_centroid("non_existent")
        assert centroid is None
    
    def test_list_families(self):
        """Test listing all families"""
        assigner = ProteinFamilyAssigner()
        
        assigner.create_family("family_001", np.random.rand(128))
        assigner.create_family("family_002", np.random.rand(128))
        
        families = assigner.list_families()
        assert len(families) == 2
        assert "family_001" in families
        assert "family_002" in families
    
    def test_get_family_stats(self):
        """Test getting family statistics"""
        assigner = ProteinFamilyAssigner()
        
        family_id = "family_001"
        assigner.create_family(family_id, np.random.rand(128))
        
        assigner.add_protein_to_family(family_id, "P12345", np.random.rand(128))
        assigner.add_protein_to_family(family_id, "P67890", np.random.rand(128))
        
        stats = assigner.get_family_stats(family_id)
        assert stats["protein_count"] == 2
        assert stats["family_id"] == family_id
    
    def test_get_family_stats_nonexistent(self):
        """Test getting stats for non-existent family"""
        assigner = ProteinFamilyAssigner()
        
        stats = assigner.get_family_stats("non_existent")
        assert stats is None
    
    def test_remove_protein_from_family(self):
        """Test removing a protein from a family"""
        assigner = ProteinFamilyAssigner()
        
        family_id = "family_001"
        protein_id = "P12345"
        
        assigner.create_family(family_id, np.random.rand(128))
        assigner.add_protein_to_family(family_id, protein_id, np.random.rand(128))
        
        assert protein_id in assigner.families[family_id]
        
        assigner.remove_protein_from_family(family_id, protein_id)
        
        assert protein_id not in assigner.families[family_id]
    
    def test_remove_family(self):
        """Test removing a family"""
        assigner = ProteinFamilyAssigner()
        
        family_id = "family_001"
        assigner.create_family(family_id, np.random.rand(128))
        
        assert family_id in assigner.families
        assert family_id in assigner.centroids
        
        assigner.remove_family(family_id)
        
        assert family_id not in assigner.families
        assert family_id not in assigner.centroids
    
    def test_update_family_centroid(self):
        """Test updating family centroid"""
        assigner = ProteinFamilyAssigner()
        
        family_id = "family_001"
        old_centroid = np.random.rand(128)
        new_centroid = np.random.rand(128)
        
        assigner.create_family(family_id, old_centroid)
        assigner.update_family_centroid(family_id, new_centroid)
        
        assert np.array_equal(assigner.centroids[family_id], new_centroid)
    
    def test_compute_similarity(self):
        """Test computing similarity between embeddings"""
        assigner = ProteinFamilyAssigner()
        
        emb1 = np.array([1.0, 0.0, 0.0])
        emb2 = np.array([0.0, 1.0, 0.0])
        emb3 = np.array([1.0, 0.0, 0.0])
        
        # Test cosine similarity
        sim1 = assigner.compute_similarity(emb1, emb2)
        sim2 = assigner.compute_similarity(emb1, emb3)
        
        assert sim1 == 0.0  # Orthogonal vectors
        assert sim2 == 1.0  # Identical vectors
    
    def test_clear_all(self):
        """Test clearing all families"""
        assigner = ProteinFamilyAssigner()
        
        assigner.create_family("family_001", np.random.rand(128))
        assigner.create_family("family_002", np.random.rand(128))
        
        assert len(assigner.families) == 2
        assert len(assigner.centroids) == 2
        
        assigner.clear_all()
        
        assert len(assigner.families) == 0
        assert len(assigner.centroids) == 0
