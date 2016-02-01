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
COPY /benchmark/GSM1480325_K562_GROseq_*.bigWig /benchmark/
COPY /benchmark/idrCode.tar.gz /benchmark/
RUN tar -zxf idrCode.tar.gz && rm idrCode.tar.gz

COPY $DIR/Makefile $DIR/
