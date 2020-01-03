#FROM brsynth/rpcache-rest
FROM brsynth/rpcache

RUN conda install -c conda-forge flask-restful

COPY rpTool.py /home/
COPY rpToolServe.py /home/
RUN rm -f /home/Dockerfile

ENTRYPOINT ["python"]
CMD ["/home/rpToolServe.py"]

# Open server port
EXPOSE 8888
