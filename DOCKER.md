# DOCKER.md - Introduction

Sattools, stvid and strf are Unix-native tools, with complex and 
specific environment requirements in order for them to operate.

The Docker environment allows you to gain access to sattools
without knowing how to install dependencies or compile source. 

It has the additional benefit of being able to use the sattools 
environment without affecting other environments on your computer.

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

Unix/Linux:

`sudo apt-get update && sudo apt-get install docker.io`

From your X11 terminal, find your computer's IP address and authorize it for X11 connections

`xauth + 192.168.1.2`

# Download and build the sattools Dockerfile
https://github.com/cbassa/sattools/blob/master/Dockerfile

## In a Terminal (Mac/Unix) or Powershell (Windows) build the Dockerfile:
`docker build  - < Dockerfile`

After a few minutes, your build should successfully finish.  Check for your image with the command:
`$ docker images`

>REPOSITORY          TAG                 IMAGE ID            CREATED             SIZE
><none>              <none>              1ad51e488e1f        12 seconds ago      1.47GB

Take note of the IMAGE ID, in this case `1ad51e488e1f` this code represents the IMAGE you just built.

To run this image in a new container, run the following command in the terminal.

```docker run -it -e DISPLAY='host.docker.internal:0 \' 
-v /tmp/.X11-unix:/tmp/.X11-unix:ro \
-v ~/Documents/satobs:/root/satellite/satobs \
--name mysattools 1ad51e488e1f \
/bin/bash
```

## Check out your new sattols Docker environment
You will be placed into the environment in the commandline of your sattools environment in the Docker container.

`root@a9031e8397c6:~# cd /root/satellite/sattools/`

The first thing you will want to do, is to specify your location. To do this, edit the file "stations.in"

`root@a9031e8397c6:~/satellite/sattools# nano data/sites.txt`

Scroll to the bottom of the file, to edit the entry for station 9999 with the name Graves.
Following the formatting for the file, enter the latitutde, longitude, elevation and name for your observing site, and save the file:

`9999 TS   47.6672 -122.0931    101    TruSat`

Next, you will want to edit the stvid configuration, which will allow you to get updated TLEs

First, register for a space-track.org account at https://www.space-track.org/auth/createAccount

Next, edit the information for your observation location, and put in the login information 
for your space-track.org account

```cd /root/satellite/stvid
root@a9031e8397c6:~/satellite/stvid# cp configuration.ini-dist configuration.ini
nano configuration.ini
```

Last, update the TLEs with

`root@a9031e8397c6:~/satellite/stvid# python3 update_tle.py`

You can see your new TLEs at:

```root@a9031e8397c6:~/satellite/stvid# cd /root/satellite/tle/
root@a9031e8397c6:~/satellite/tle# ls -al
```

Here a few things that you can do now that your environment is set up for your location with TLEs

Get sunset / sunrise for your site (in UTC):

`root@a9031e8397c6:~/satellite/tle# allnight`

> 2019-11-27T00:56:27 2019-11-27T14:55:44

Run passes for your site for 3 hours, above 45 degrees, starting at sunset (from above):
`pass -t 2019-11-27T00:56:27 -l 10800  -A 45`

Run skymap to interactively view passes for your site:
`skymap S`

This should open a graphical display in an X-window on your system.
Use > to advance by 1 minute

Type `exit` to leave the docker environment.

# Returning to your docker environment.

You can check to see if your container is still running by typing the following in the terminal:
```(MVP) MBP-140320:~ chris$ docker ps -a
CONTAINER ID        IMAGE               COMMAND                  CREATED             STATUS                  PORTS               NAMES
a9031e8397c6        d82642132ae1        "/bin/sh -c /bin/bas…"   6 days ago          Up 6 days                                   st-b3
6fe73f3d652e        d82642132ae1        "/bin/sh -c /bin/bas…"   6 days ago          Exited (0) 6 days ago                       st-b2
```

You can see multiple containers from the image you built.  In this case, the container a9031e8397c6 is still running.

To reattach to this container, execute:

`docker attach a9031e8397c6`

You might need to press `return` to see the prompt.

If you don't have a container running (like container 6fe73f3d652e), you can restart it with:

`docker start 6fe73f3d652e`

If you don't have any active containers, you can start over with the `docker run` command above.

## Things to watch out for
1. When you build/start/stop/exit containers, its possible for multiple (1.5 gigabyte) containers to pile up.  Review the commands for `docker rmi` and `docker rm` for how to remove thiese.
1. Any changes you make inside the docker container are not automatively saved to the image. If you would like to save your changes, use the `docker commit` command to make an updated version of the images with your changes.

