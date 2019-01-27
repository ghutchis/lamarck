"""
Unit and regression test for the lamarck package.
"""

# Import package, test suite, and other packages as needed
import lamarck
import pytest
import sys

def test_lamarck_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "lamarck" in sys.modules
