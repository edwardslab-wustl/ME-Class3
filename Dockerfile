FROM centos:7.7.1908

LABEL maintainer="Manoj Singh <manoj@wustl.edu>"

RUN yum install -y git python3 numpy scipy scikit-learn python-matplotlib ipython python-pandas sympy python-nose atlas-devel && \
    pip3 install Pillow==6.2.2 && \
    cd tmp && pip3 install git+https://github.com/edwardslab-wustl/ME-Class3.git && \
    yum clean all 
#    rm -fr /tmp/ME-Class3
WORKDIR /opt
RUN git clone https://github.com/edwardslab-wustl/ME-Class3.git

