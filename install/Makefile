IPS_DIR = /global/common/software/atom/cori

PREFIX = ${IPS_DIR}/cesol/ips_fastran

RSYNC = rsync -avz
OPTIONS = --exclude "__pycache__"
DIRS = src

VERSION_FILE = $(IPS_DIR)/cesol/VERSION.IPS_FASTRAN

ifndef LIB_DIR 
   LIB_DIR = ${IPS_DIR}/cesol/share
endif
ifndef CONF_DIR 
   CONF_DIR = ${IPS_DIR}/cesol/conf
endif

poo:
	@echo 'nothing done'

all: dir install version

dir:
	@echo 'make build dir:' $(PREFIX)
	@test -d $(PREFIX) || mkdir -p $(PREFIX)
	@test -d $(PREFIX)/src || mkdir -p $(PREFIX)/src
	@test -d $(PREFIX)/bin || mkdir -p $(PREFIX)/bin
	@test -d $(PREFIX)/doc || mkdir -p $(PREFIX)/doc
	@test -d $(LIB_DIR) || mkdir -p $(LIB_DIR)
	@test -d $(CONF_DIR) || mkdir -p $(CONF_DIR)
	@echo $(DIRS)

install: ${DIRS}
	@echo ${PREFIX}
	${RSYNC} ${OPTIONS} --delete src/ ${PREFIX}/src/
	${RSYNC} ${OPTIONS} lib/ ${LIB_DIR}/
ifneq ($(LIB_DIR), $(PREFIX)/lib)
	rm -f $(PREFIX)/lib
	ln -s ${LIB_DIR} $(PREFIX)/lib 
endif
	${RSYNC} ${OPTIONS} conf/ ${CONF_DIR}/

version:
	git log --format="%H" -n 1 > $(VERSION_FILE)
	date >> $(VERSION_FILE)

