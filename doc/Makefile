

examples = $(wildcard example-*)
exoutput = $(foreach f,$(examples),$(f)/_output)

.PHONY: html
#sphinx: $(exoutput)
html: 
	sphinx-build -b html . ./users

# to make output needed for plot examples:
$(exoutput): 
	@echo 'Making ' $@
	$(MAKE) --directory=$(subst /_output,,$@) .rst


