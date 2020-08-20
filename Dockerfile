FROM dask-base:latest
WORKDIR /usr/src/app

COPY ./requirements.txt /usr/src/app
RUN pip install -r requirements.txt

COPY ./src /usr/src/app
RUN python setup_cython.py build_ext --inplace