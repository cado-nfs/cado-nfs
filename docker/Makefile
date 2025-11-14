SHELL:=bash

all: \
	docker-image-factoring-full	\
	docker-image-factoring-client	\
	docker-image-factoring-server	\
	docker-image-factoring-linalg

# The discretelog images are not ready yet, we should think about the
# best way to do it.

NOW:=$(shell date +%Y%m%d%H%M)

# REG=registry.gitlab.inria.fr/cado-nfs/cado-nfs
REG:=$$CI_REGISTRY_IMAGE

docker-image-%: docker/base docker/%
	./tests/provide-wdir.sh --env WDIR bash -xc 'git archive --format=tar.gz HEAD > $$WDIR/cado.tar.gz ; cp -p docker/local.sh docker/dispatch-files.sh $$WDIR/ ; cat $^ > $$WDIR/Dockerfile ; docker build -t $(REG)/$* $$WDIR'
	! [ "$$PUSH" ] || ! [ "$$CI_COMMIT_BRANCH" ] || [ "$$CI_COMMIT_BRANCH" != master ] || (docker tag $*:latest $(REG)/$*:$(NOW) ; docker push $(REG)/$*:$(NOW) ; docker push $(REG)/$*)

