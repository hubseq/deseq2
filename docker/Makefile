NAME=deseq2
TAG=00.00.03

all: build push

build:
	docker build -t $(NAME):$(TAG) -f Dockerfile ../
	docker tag $(NAME):$(TAG) $(REGISTRY):$(TAG)

push:
	docker push $(REGISTRY):$(TAG)
