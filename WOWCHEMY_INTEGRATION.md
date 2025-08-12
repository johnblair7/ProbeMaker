# Integrating ProbeMaker into Your Wowchemy Website

This guide will help you integrate the ProbeMaker tool into your existing Wowchemy (Hugo Academic) website.

## 🎯 **Integration Options**

### **Option 1: Simple HTML Embed (Easiest)**
Embed the standalone HTML file directly into a Wowchemy page.

### **Option 2: Full Backend Integration (Most Functional)**
Deploy the Flask backend and connect it to your website.

### **Option 3: Hybrid Approach (Recommended)**
Use the standalone HTML for the interface, with a separate backend service.

---

## **Option 1: Simple HTML Embed**

### **Step 1: Add to Your Wowchemy Content**
1. Copy `probe_maker_standalone.html` to your Wowchemy content directory
2. Create a new page in your Wowchemy site
3. Embed the HTML using Hugo's shortcodes or raw HTML

### **Step 2: Create a Wowchemy Page**
Create a new file in your content directory (e.g., `content/tools/probe-maker.md`):

```markdown
---
title: "ProbeMaker - mRNA Probe Generator"
subtitle: "Generate complementary sequences to mRNA for your genes"
summary: "A web-based tool for designing DNA probes complementary to mRNA sequences"
authors: [your-name]
tags: [tools, bioinformatics, molecular-biology]
categories: [tools]
date: 2024-01-01
lastmod: 2024-01-01
featured: false
draft: false

# Featured image
image:
  caption: "ProbeMaker Interface"
  focal_point: ""
  preview_only: false

# Projects (optional).
projects: []
---

## About ProbeMaker

ProbeMaker is a tool for generating complementary DNA sequences to mRNA. It automatically fetches gene sequences from NCBI and designs optimal probes with specific constraints.

{{< rawhtml >}}
{{< read_file path="probe_maker_standalone.html" >}}
{{< /rawhtml >}}
```

### **Step 3: Customize the Styling**
The standalone HTML includes its own CSS that should work well with most Wowchemy themes. You can customize colors and fonts to match your site's theme.

---

## **Option 2: Full Backend Integration**

### **Step 1: Deploy the Flask Backend**
Follow the [DEPLOYMENT.md](DEPLOYMENT.md) guide to deploy the Flask app to your hosting provider.

### **Step 2: Update the HTML to Point to Your Backend**
Modify the JavaScript in `probe_maker_standalone.html` to point to your deployed backend:

```javascript
// Replace this line in the fetch call:
const response = await fetch('/generate', {
    method: 'POST',
    body: formData
});

// With your backend URL:
const response = await fetch('https://your-backend-domain.com/generate', {
    method: 'POST',
    body: formData
});
```

### **Step 3: Handle CORS (if needed)**
If your backend is on a different domain, you may need to enable CORS in your Flask app:

```python
from flask_cors import CORS

app = Flask(__name__)
CORS(app)  # Enable CORS for all routes
```

---

## **Option 3: Hybrid Approach (Recommended)**

### **Why This Approach?**
- **Frontend**: Beautiful interface integrated into your Wowchemy site
- **Backend**: Separate service that can be deployed anywhere
- **Flexibility**: Easy to maintain and update independently

### **Implementation:**
1. Use the standalone HTML in your Wowchemy site
2. Deploy the Flask backend to a service like:
   - **Heroku** (free tier available)
   - **Railway** (simple deployment)
   - **Render** (free tier available)
   - **Your own server** (if you have one)

---

## **🎨 Customization for Wowchemy**

### **Theme Integration**
The standalone HTML uses neutral colors that should work with most Wowchemy themes. You can customize:

```css
/* Update these CSS variables to match your theme */
:root {
    --primary-color: #3498db;      /* Your theme's primary color */
    --secondary-color: #2980b9;    /* Your theme's secondary color */
    --background-color: #f8f9fa;   /* Your theme's background */
    --text-color: #2c3e50;         /* Your theme's text color */
}
```

### **Font Integration**
The HTML uses system fonts by default. To use your Wowchemy theme's fonts:

```css
body {
    font-family: var(--font-family-base), -apple-system, BlinkMacSystemFont, sans-serif;
}
```

---

## **📱 Mobile Responsiveness**

The standalone HTML is fully responsive and will work well on mobile devices. It automatically adjusts:
- Font sizes
- Layout spacing
- Button sizes
- Form elements

---

## **🔧 Advanced Customization**

### **Adding Analytics**
If you want to track usage, add Google Analytics or similar:

```html
<!-- Add this in the <head> section -->
<script async src="https://www.googletagmanager.com/gtag/js?id=GA_MEASUREMENT_ID"></script>
<script>
  window.dataLayer = window.dataLayer || [];
  function gtag(){dataLayer.push(arguments);}
  gtag('js', new Date());
  gtag('config', 'GA_MEASUREMENT_ID');
</script>
```

### **Adding Social Sharing**
Add social media sharing buttons below the tool:

```html
<div class="social-sharing">
    <button onclick="shareOnTwitter()">Share on Twitter</button>
    <button onclick="shareOnLinkedIn()">Share on LinkedIn</button>
</div>
```

---

## **🚀 Deployment Checklist**

- [ ] Choose your integration approach
- [ ] Test the standalone HTML locally
- [ ] Deploy backend (if using full integration)
- [ ] Update API endpoints in HTML
- [ ] Test on your Wowchemy site
- [ ] Customize styling to match your theme
- [ ] Add analytics (optional)
- [ ] Test on mobile devices

---

## **❓ Troubleshooting**

### **Common Issues:**

1. **Styling Conflicts**: The standalone HTML includes its own CSS that should override most conflicts
2. **JavaScript Errors**: Check browser console for any JavaScript errors
3. **API Connection**: Ensure your backend URL is correct and accessible
4. **Mobile Issues**: Test on various devices and screen sizes

### **Getting Help:**
- Check the browser console for error messages
- Verify your backend is running and accessible
- Test the standalone HTML file independently first

---

## **🎉 Success!**

Once integrated, your Wowchemy site will have a professional, functional tool that:
- Looks great on all devices
- Integrates seamlessly with your site's design
- Provides real value to your visitors
- Demonstrates your technical expertise

The tool will be accessible to anyone visiting your website, making your research tools available to the broader scientific community!
