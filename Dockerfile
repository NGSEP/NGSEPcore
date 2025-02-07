FROM openjdk:11.0.8-slim

WORKDIR /app

RUN mkdir build

COPY . build/

RUN apt-get update
RUN apt-get install -y make

RUN cd build && make all
RUN cp ./build/*.jar /app/
RUN rm -rf ./build

CMD [ "sh" ]
