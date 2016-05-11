FROM nanakiksc/lab_docker

# Install other software. Zerone is already in the base image.
RUN echo 'source("http://bioconductor.org/biocLite.R")' > install.R && \
    echo 'biocLite("BayesPeak")' >> install.R && \
    mkdir /BayesPeak && mv install.R /BayesPeak && \
    R -f /BayesPeak/install.R
RUN git clone https://github.com/mahmoudibrahim/JAMM.git && \
    echo 'install.packages(c("signal", "mclust"), repos="http://cran.us.r-project.org")' > /JAMM/dependencies.R && \
    R -f /JAMM/dependencies.R
# MACS depends on NumPy.
RUN apt-get update && apt-get install -y \
    python-numpy=1:1.8.2-0* && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    wget https://bootstrap.pypa.io/get-pip.py -O - | python && \
    pip install MACS2
ENV PATH=$PATH:/home/bin \
    PYTHONPATH=$PYTHONPATH:/home/lib/python2.7/site-packages
RUN wget https://github.com/bedops/bedops/releases/download/v2.4.14/bedops_linux_x86_64-v2.4.14.tar.bz2 && \
    tar xjf bedops_linux_x86_64-v2.4.14.tar.bz2 bin/sam2bed && \
    rm bedops_linux_x86_64-v2.4.14.tar.bz2

ENV DIR /discretize
RUN mkdir $DIR
COPY $DIR/datasets.md5 $DIR/
COPY $DIR/ignore_random_chrom.py $DIR/
COPY /BayesPeak/BayesPeak2bed.py /BayesPeak/
COPY /BayesPeak/BayesPeak.R /BayesPeak/

RUN mkdir /benchmark
COPY /benchmark/Makefile /benchmark/
COPY /benchmark/ctcf_motif_coverage.R /benchmark/
COPY /benchmark/ctcf_motifs.bed /benchmark/
COPY /benchmark/tss_coverage.R /benchmark/
COPY /benchmark/hg19_TSS.bed /benchmark/
COPY /benchmark/performance.R /benchmark/
COPY /benchmark/GSM1480325_K562_GROseq_minus.bigWig /benchmark/
COPY /benchmark/GSM1480325_K562_GROseq_plus.bigWig /benchmark/
COPY /benchmark/H1ESC_files_path.txt /benchmark/
COPY /benchmark/H1ESC_files_pairing.txt /benchmark/
COPY /benchmark/idrCode.tar.gz /benchmark/
RUN tar xzf /benchmark/idrCode.tar.gz && rm /benchmark/idrCode.tar.gz

COPY $DIR/Makefile $DIR/
