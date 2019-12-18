FROM brsynth/rpbase

#RUN apt-get update && apt-get -y upgrade
RUN apt-get install -y lzma

COPY rpTool.py /home/
COPY rpToolServe.py /home/
