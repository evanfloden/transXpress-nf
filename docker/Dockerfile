FROM continuumio/miniconda
MAINTAINER Chris Wyatt <chris-wyatt@hotmail.com>

RUN apt-get -y install ttf-dejavu

COPY conda.yml .
RUN \
   conda env update -n root -f conda.yml \
&& conda clean -a
