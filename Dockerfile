# Use the official Python 3.11 image as a parent image
FROM python:3.11-rc-slim

# Set the working directory in the container
WORKDIR /usr/src/app

# Install poetry
RUN pip install --no-cache-dir poetry

# Copy only the pyproject.toml and (if exists) poetry.lock files into the working directory
COPY pyproject.toml poetry.lock* /usr/src/app/

# Disable virtualenv creation by poetry and install dependencies
# Virtualenvs are not useful for Docker since the container itself is already isolated.
RUN poetry config virtualenvs.create false \
  && poetry install --no-interaction --no-ansi

# Copy the rest of your application code
COPY . .

# Command to run on container start
CMD [ "python", "./filter_vcf_and_cadd_scores/main.py" ]
