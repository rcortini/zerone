FROM nanakiksc/lab_docker

# Install other software. Zerone is already in the base image.
RUN mkdir /BayesPeak && \
    echo \
    $'source("https://bioconductor.org/biocLite.R")\nbiocLite("BayesPeak")' \
    > /BayesPeak/install.R && \
    R -f /BayesPeak/install.R
RUN git clone https://github.com/mahmoudibrahim/JAMM.git && \
    echo 'install.packages(c("signal", "mclust"))' > /JAMM/dependencies.R && \
    R -f /JAMM/dependencies.R
# MACS depends on NumPy.
RUN apt-get update && apt-get install -y \
    python-numpy=1:1.8.2-0* && \
    apt-get clean && rm -rf /var/lib/apt/lists/* && \
    git clone https://github.com/taoliu/MACS.git && \
    python MACS/setup.py install
ENV PATH=$PATH:/home/bin \
    PYTHONPATH=$PYTHONPATH:/home/lib/python2.7/site-packages

ENV DIR /replay
RUN mkdir $DIR
COPY datasets.md5 $DIR/
COPY ignore_random_chrom.py $DIR/
COPY Makefile $DIR/
COPY BayesPeak2bed.py /BayesPeak
