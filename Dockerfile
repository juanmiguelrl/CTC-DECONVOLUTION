FROM python:3.10
COPY requirements.txt requirements.txt
RUN python -m pip install -r requirements.txt
COPY src .
ENTRYPOINT ["python", "main.py"]