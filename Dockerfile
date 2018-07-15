FROM python:3.6

RUN apt update && apt install -y \
    python3-dev \
    python3-setuptools \
    python3-sphinx

RUN apt update && apt install -y \
    texlive-generic-recommended \
    texlive-generic-extra \
    texlive-fonts-recommended \
    texlive-latex-recommended \
    texlive-latex-extra \
    texlive-extra-utils
