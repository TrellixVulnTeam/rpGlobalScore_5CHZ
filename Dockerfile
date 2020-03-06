FROM brsynth/rpbase:dev

RUN pip install scipy

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY tool_rpGlobalScore.py /home/
