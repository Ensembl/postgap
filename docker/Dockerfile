FROM ubuntu

RUN apt-get update && apt-get -y install build-essential curl git vim wget unzip zlib1g-dev python2.7 python-pip autoconf dh-autoreconf pkg-config libatlas-base-dev

# these python packages also feature in install_dependencies.sh
# but we install them here as root to avoid sudo-ing or path issues
RUN pip install pybedtools requests pandas flask cherrypy

# create a postgap user
RUN useradd -r -m -U -d /home/postgap -s /bin/bash -c "POSTGAP User" -p '' postgap
RUN usermod -a -G sudo postgap
USER postgap
ENV HOME /home/postgap

# create src dir
RUN mkdir -p $HOME/src
WORKDIR $HOME/src

# fetch and set up postgap
RUN git clone https://github.com/Ensembl/postgap.git

WORKDIR postgap

USER postgap
RUN sh install_dependencies.sh

# set up env
ENV PYTHONPATH $HOME/src/postgap/lib/
ENV PATH $PATH:$HOME/src/postgap/bin/
