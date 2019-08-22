PFPR=PfPR/Raster\ Data/PfPR_median/PfPR_median_Global_admin0_2017.tif
HRSL=hrsl/hrsl_uga_settlement.tif

all:  $(HRSL) $(PFPR)

hrsl_uga_v1.zip:
	curl -o hrsl_uga_v1.zip https://www.ciesin.columbia.edu/repository/hrsl/hrsl_uga_v1.zip

$(HRSL): hrsl_uga_v1.zip
	mkdir -p hrsl
	(cd hrsl && unzip -n ../hrsl_uga_v1.zip)
	touch $(HRSL)

PfPR.zip:
	curl -o PfPR.zip https://map.ox.ac.uk/wp-content/uploads/lancet-2019/pf/PfPR.zip

$(PFPR): PfPR.zip
	unzip -n PfPR.zip
	touch $(PFPR)

clean:
	rm -rf hrsl PfPR *.zip
