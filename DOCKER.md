# DOCKER.md - Introduction

Sattools, stvid and strf are Unix-native tools, with complex and 
specific environment requirements in order for them to operate.

The Docker environment allows you to gain access to sattools
without knowing how install dependencies, compile or use the 
sattools environment without affecting other environments on 
your computer.

The instructions below will get you running sattools in about 
15 minutes.

# Setup your Docker / X-Windows environment

Mac/PC **Docker Desktop**
- Create an account and download Docker Desktop:
  https://www.docker.com/products/docker-desktop

Windows:
    Install **Xming X Server for Windows** - https://sourceforge.net/projects/xming/
Mac:
    Install **XQuartz** - https://www.xquartz.org/ (can also install via brew)
    Log out/back in as directed by the installer

Unix/Linux
`sudo apt-get install docker.io`

From your X11 terminal, find your computer's IP address and authorize it for X11 connections
`xauth + 192.168.1.2`
