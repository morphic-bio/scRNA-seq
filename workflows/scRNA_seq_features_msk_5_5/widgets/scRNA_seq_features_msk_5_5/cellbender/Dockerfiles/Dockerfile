FROM cumulusprod/cellbender:0.3.0
RUN apt-get update && apt-get install -y --no-install-recommends git  && rm -rf /var/lib/apt/lists/*
RUN git clone https://github.com/broadinstitute/CellBender.git
RUN pip install -e CellBender
COPY remove_noise.sh /usr/local/bin/.
#patch for copying jupyter notebook
COPY report.py CellBender/cellbender/remove_background/report.py
COPY addCounts.py /usr/local/bin/.

