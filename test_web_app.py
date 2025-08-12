#!/usr/bin/env python3
"""
Simple test script for the ProbeMaker web application.
Run this to test the web endpoints without the full web interface.
"""

import requests
import time

def test_web_app():
    """Test the web app endpoints."""
    base_url = "http://localhost:5000"
    
    print("🧪 Testing ProbeMaker Web App...")
    
    # Test health endpoint
    try:
        response = requests.get(f"{base_url}/health")
        if response.status_code == 200:
            print("✅ Health check: PASSED")
        else:
            print(f"❌ Health check: FAILED - Status {response.status_code}")
    except requests.exceptions.ConnectionError:
        print("❌ Health check: FAILED - Cannot connect to server")
        print("   Make sure the web app is running with: python app.py")
        return
    
    # Test probe generation with a few genes
    test_genes = "BRCA1\nTP53\nEGFR"
    
    print("\n🧬 Testing probe generation...")
    try:
        response = requests.post(f"{base_url}/generate", data={'gene_names': test_genes})
        
        if response.status_code == 200:
            result = response.json()
            if result.get('success'):
                print("✅ Probe generation: PASSED")
                print(f"   Generated {result.get('probe_count', 0)} probe pairs")
                print(f"   For {result.get('gene_count', 0)} genes")
                
                # Test file download
                file_path = result.get('file_path', '').split('/')[-1]
                if file_path:
                    print(f"\n📥 Testing file download...")
                    download_response = requests.get(f"{base_url}/download/{file_path}")
                    if download_response.status_code == 200:
                        print("✅ File download: PASSED")
                        print(f"   File size: {len(download_response.content)} bytes")
                    else:
                        print(f"❌ File download: FAILED - Status {download_response.status_code}")
            else:
                print(f"❌ Probe generation: FAILED - {result.get('error', 'Unknown error')}")
        else:
            print(f"❌ Probe generation: FAILED - Status {response.status_code}")
            print(f"   Response: {response.text}")
            
    except Exception as e:
        print(f"❌ Probe generation: FAILED - Exception: {e}")
    
    print("\n🎉 Web app testing complete!")
    print("\nTo use the web interface:")
    print(f"1. Open your browser to: {base_url}")
    print("2. Enter gene names in the text area")
    print("3. Click 'Generate Probe Sequences'")
    print("4. Download the results file")

if __name__ == "__main__":
    test_web_app()
