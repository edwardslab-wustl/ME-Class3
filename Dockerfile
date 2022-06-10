FROM centos:7.7.1908

LABEL maintainer="John Edwards <jredwards@wustl.edu>"

RUN yum install -y git python3 numpy scipy scikit-learn python-matplotlib ipython python-pandas sympy python-nose atlas-devel && \
    yum clean all 

RUN pip3 install Pillow==6.2.2
RUN cd tmp && pip3 install git+https://github.com/edwardslab-wustl/ME-Class3.git

WORKDIR /opt
RUN git clone https://github.com/edwardslab-wustl/ME-Class3.git

