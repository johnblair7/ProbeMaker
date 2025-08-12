# ProbeMaker Web App Deployment Guide

This guide will help you deploy the ProbeMaker web application on your website.

## Prerequisites

- Python 3.7+
- pip (Python package manager)
- A web server (Apache/Nginx) or cloud hosting service

## Option 1: Simple Development Server (Testing)

For testing or development purposes:

```bash
# Install dependencies
pip install -r requirements.txt

# Run the development server
python app.py
```

The app will be available at `http://localhost:5000`

## Option 2: Production Deployment with Gunicorn

### 1. Install Dependencies

```bash
pip install -r requirements.txt
pip install gunicorn
```

### 2. Run with Gunicorn

```bash
# Basic gunicorn command
gunicorn -c gunicorn.conf.py app:app

# Or with custom settings
gunicorn --bind 0.0.0.0:8000 --workers 4 --timeout 300 app:app
```

### 3. Using Systemd (Linux)

Create a systemd service file `/etc/systemd/system/probe_maker.service`:

```ini
[Unit]
Description=ProbeMaker Web Application
After=network.target

[Service]
User=www-data
Group=www-data
WorkingDirectory=/path/to/your/probe_maker
Environment="PATH=/path/to/your/venv/bin"
ExecStart=/path/to/your/venv/bin/gunicorn -c gunicorn.conf.py app:app
Restart=always

[Install]
WantedBy=multi-user.target
```

Enable and start the service:

```bash
sudo systemctl enable probe_maker
sudo systemctl start probe_maker
sudo systemctl status probe_maker
```

## Option 3: Cloud Deployment

### Heroku

1. Create a `Procfile`:
```
web: gunicorn -c gunicorn.conf.py app:app
```

2. Deploy:
```bash
heroku create your-app-name
git push heroku main
```

### Docker

Create a `Dockerfile`:

```dockerfile
FROM python:3.9-slim

WORKDIR /app
COPY requirements.txt .
RUN pip install -r requirements.txt

COPY . .
EXPOSE 8000

CMD ["gunicorn", "-c", "gunicorn.conf.py", "app:app"]
```

Build and run:
```bash
docker build -t probe_maker .
docker run -p 8000:8000 probe_maker
```

## Option 4: Shared Hosting (cPanel, etc.)

1. Upload all files to your hosting directory
2. Set up a Python app in your hosting control panel
3. Point it to `app.py`
4. Install dependencies via pip or requirements.txt

## Environment Variables

Set these environment variables for production:

```bash
export SECRET_KEY="your-secret-key-here"
export FLASK_ENV="production"
```

## Reverse Proxy with Nginx

If using Nginx as a reverse proxy, add this to your Nginx config:

```nginx
server {
    listen 80;
    server_name yourdomain.com;

    location / {
        proxy_pass http://127.0.0.1:8000;
        proxy_set_header Host $host;
        proxy_set_header X-Real-IP $remote_addr;
        proxy_set_header X-Forwarded-For $proxy_add_x_forwarded_for;
        proxy_set_header X-Forwarded-Proto $scheme;
        
        # Increase timeout for NCBI requests
        proxy_read_timeout 300;
        proxy_connect_timeout 300;
        proxy_send_timeout 300;
    }
}
```

## Security Considerations

1. **Rate Limiting**: Consider adding rate limiting to prevent abuse
2. **File Cleanup**: Temporary files are automatically cleaned up
3. **Input Validation**: Gene names are validated and limited to 100 per request
4. **HTTPS**: Use HTTPS in production for security

## Monitoring

The app includes a health check endpoint at `/health` for monitoring:

```bash
curl https://yourdomain.com/health
```

## Troubleshooting

### Common Issues

1. **NCBI Rate Limiting**: The app handles this automatically with delays
2. **Memory Issues**: Reduce worker count in gunicorn.conf.py
3. **Timeout Errors**: Increase timeout values for slow NCBI responses

### Logs

Check logs for errors:
```bash
# If using systemd
sudo journalctl -u probe_maker -f

# If using gunicorn directly
gunicorn --log-level debug app:app
```

## Performance Tips

1. **Worker Count**: Adjust based on your server's CPU cores
2. **Caching**: Consider adding Redis for caching gene sequences
3. **CDN**: Use a CDN for static assets if serving many users

## Support

For issues or questions, check the main README.md or create an issue in the repository.
