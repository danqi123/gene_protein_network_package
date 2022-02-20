FROM python:3.9.9

WORKDIR .
RUN pip install --upgrade pip
COPY . .
RUN pip install -e plab2_package

LABEL "author"="Danqi" "email"="wdanqi@live.com"
LABEL build_date="2021-12-19"

#FLASK_PORT 5005
ENV FLASK_DEBUG=true FLASK_APP=run.py FLASK_ENV=development
ENV FLASK_PORT=5005
EXPOSE $FLASK_PORT

ENTRYPOINT [ "python", "flask_template/run.py" ]




