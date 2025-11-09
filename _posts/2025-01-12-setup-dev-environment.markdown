---
layout: post
title: "Setting up the perfect development environment"
date: 2025-01-12 14:30:00 +0100
categories: tutorial devops
excerpt: "How to configure a modern development environment with Docker, VS Code and essential tools."
---

# Setting up the perfect development environment

A well-configured development environment can make a huge difference in daily productivity.

## Essential tools

### Editor
- **VS Code** with specific extensions
- **Vim** for quick editing on servers
- **JetBrains IDEs** for complex projects

### Terminal
```bash
# Oh My Zsh with useful plugins
sh -c "$(curl -fsSL https://raw.github.com/ohmyzsh/ohmyzsh/master/tools/install.sh)"
```

### Containerization
Docker is now essential:
```dockerfile
FROM node:18-alpine
WORKDIR /app
COPY package*.json ./
RUN npm ci --only=production
COPY . .
EXPOSE 3000
CMD ["npm", "start"]
```

## Conclusions

A well-thought setup pays off over time.