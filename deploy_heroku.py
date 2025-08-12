#!/usr/bin/env python3
"""
Simple Heroku deployment script for ProbeMaker web app.
Run this to quickly deploy to Heroku.
"""

import os
import subprocess
import sys

def run_command(command, description):
    """Run a shell command and handle errors."""
    print(f"🔄 {description}...")
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        print(f"✅ {description} completed successfully")
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed:")
        print(f"   Error: {e.stderr}")
        return None

def create_heroku_files():
    """Create necessary Heroku deployment files."""
    print("📝 Creating Heroku deployment files...")
    
    # Create Procfile
    with open('Procfile', 'w') as f:
        f.write('web: gunicorn -c gunicorn.conf.py app:app\n')
    
    # Create runtime.txt
    with open('runtime.txt', 'w') as f:
        f.write('python-3.9.18\n')
    
    # Create app.json for Heroku
    app_json = '''{
  "name": "probe-maker",
  "description": "Generate complementary sequences to mRNA for given genes",
  "repository": "https://github.com/yourusername/probe-maker",
  "logo": "https://node-js-sample.herokuapp.com/node.png",
  "keywords": ["python", "flask", "bioinformatics", "molecular-biology"],
  "env": {
    "SECRET_KEY": {
      "description": "A secret key for the Flask application",
      "generator": "secret"
    }
  },
  "buildpacks": [
    {
      "url": "heroku/python"
    }
  ]
}'''
    
    with open('app.json', 'w') as f:
        f.write(app_json)
    
    print("✅ Heroku files created")

def deploy_to_heroku():
    """Deploy the app to Heroku."""
    print("🚀 Starting Heroku deployment...")
    
    # Check if Heroku CLI is installed
    if not run_command("heroku --version", "Checking Heroku CLI"):
        print("❌ Heroku CLI not found. Please install it first:")
        print("   https://devcenter.heroku.com/articles/heroku-cli")
        return False
    
    # Check if git repository exists
    if not os.path.exists('.git'):
        print("❌ Not a git repository. Please run 'git init' first.")
        return False
    
    # Create Heroku app
    app_name = input("Enter your Heroku app name (or press Enter for auto-generated): ").strip()
    
    if app_name:
        create_result = run_command(f"heroku create {app_name}", f"Creating Heroku app '{app_name}'")
    else:
        create_result = run_command("heroku create", "Creating Heroku app")
    
    if not create_result:
        return False
    
    # Get the app URL from the output
    app_url = None
    for line in create_result.split('\n'):
        if 'https://' in line and 'herokuapp.com' in line:
            app_url = line.strip()
            break
    
    if not app_url:
        print("❌ Could not determine Heroku app URL")
        return False
    
    print(f"🌐 Your app will be available at: {app_url}")
    
    # Add and commit files
    run_command("git add .", "Adding files to git")
    run_command("git commit -m 'Deploy to Heroku'", "Committing changes")
    
    # Push to Heroku
    if not run_command("git push heroku main", "Pushing to Heroku"):
        return False
    
    # Open the app
    run_command("heroku open", "Opening app in browser")
    
    print(f"\n🎉 Deployment complete!")
    print(f"🌐 Your app is live at: {app_url}")
    print(f"📊 Monitor logs with: heroku logs --tail")
    print(f"🔧 Open Heroku dashboard: heroku dashboard")
    
    return True

def main():
    """Main deployment function."""
    print("🧬 ProbeMaker Heroku Deployment")
    print("=" * 40)
    
    # Check prerequisites
    if not os.path.exists('app.py'):
        print("❌ app.py not found. Please run this script from the ProbeMaker directory.")
        return
    
    if not os.path.exists('requirements.txt'):
        print("❌ requirements.txt not found. Please run this script from the ProbeMaker directory.")
        return
    
    # Create Heroku files
    create_heroku_files()
    
    # Deploy
    if deploy_to_heroku():
        print("\n🎊 Congratulations! Your ProbeMaker app is now live on Heroku!")
        print("\nNext steps:")
        print("1. Update the API endpoint in your Wowchemy HTML file")
        print("2. Test the app functionality")
        print("3. Customize the styling if needed")
    else:
        print("\n❌ Deployment failed. Please check the error messages above.")

if __name__ == "__main__":
    main()
