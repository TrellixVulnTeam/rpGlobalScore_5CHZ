FROM brsynth/rpbase

RUN pip install scipy

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY galaxy/code/tool_rpGlobalScore.py /home/
