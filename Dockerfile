FROM kbase/sdkbase2:python
MAINTAINER KBase Developer
# -----------------------------------------

# Insert apt-get instructions here to install
# any required dependencies for your module.

# RUN apt-get update
RUN apt-get update -qq && \
    apt-get install -yq --no-install-recommends \
                                                git \
                                                less \
                                                libdatetime-perl \
                                                libxml-simple-perl \
                                                libdigest-md5-perl \
                                                bioperl && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# clone prokka
WORKDIR /kb
RUN git clone -b "v1.12" https://github.com/tseemann/prokka && \
    prokka/bin/prokka --setupdb

# set links to /usr/bin
ENV PATH $PATH:/kb/prokka/bin

# Update tbl2asn to recent version (25.3):
#RUN mv /kb/prokka/binaries/linux/tbl2asn /kb/prokka/binaries/linux/tbl2asn.orig && \
#    wget ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz && \
#    gunzip ./linux.tbl2asn.gz && \
#    chmod 777 ./linux.tbl2asn && \
#    mv ./linux.tbl2asn /kb/prokka/binaries/linux/tbl2asn

WORKDIR /kb
RUN git clone -b "0.8" https://github.com/tseemann/barrnap
ENV PATH $PATH:/kb/barrnap/bin

RUN git clone -b "bcbio-gff-v0.6.4" https://github.com/chapmanb/bcbb && \
    cd bcbb/gff && \
    python setup.py build && \
    python setup.py install

# -----------------------------------------

COPY ./ /kb/module
RUN mkdir -p /kb/module/work
RUN chmod -R a+rw /kb/module

WORKDIR /kb/module

RUN make all

ENTRYPOINT [ "./scripts/entrypoint.sh" ]

CMD [ ]
