FROM brsynth/rpbase:dev

RUN pip install scipy

COPY rpTool.py /home/
COPY rpToolServe.py /home/
