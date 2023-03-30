FROM continuumio/miniconda3:22.11.1

LABEL maintainer="John Edwards <jredwards@wustl.edu>"

# Install python3 and other needed rpms
RUN apt-get update \
    && apt-get install -y --no-install-recommends git \
    && apt-get clean
#RUN yum install -y git python3 python3-devel atlas-devel zlib-devel \
#        libjpeg-turbo-devel gcc libtiff-devel freetype-devel tcl-devel \
#        libwebp-devel \
#        && yum clean all 

# dev packages (remove later)
#RUN yum install -y epel-release
#RUN yum update -y
#RUN yum install -y htop \
#    && yum clean all

#Install python packages
#RUN pip3 install numpy scipy pandas matplotlib dataclasses scikit-learn
#RUN pip3 install numpy scipy pandas matplotlib dataclasses
#RUN pip3 install seaborn
WORKDIR /usr/local
COPY requirements.txt /usr/local/requirements.txt
RUN pip3 install -r requirements.txt
RUN rm requirements.txt

#Install ME-Class3
WORKDIR /usr/local
ADD https://api.github.com/repos/edwardslab-wustl/ME-Class3/git/refs/heads/main version.json
RUN mkdir ME-Class3
COPY . /usr/local/ME-Class3

WORKDIR /usr/local/ME-Class3
#RUN pip3 install --use-feature=in-tree-build -e .
#RUN pip3 install -r requirements.txt
RUN pip3 install -e .
RUN chmod -R 777 /usr/local/ME-Class3
#RUN cd tmp && pip3 install git+https://github.com/edwardslab-wustl/ME-Class3.git
#RUN git clone https://github.com/edwardslab-wustl/ME-Class3.git

