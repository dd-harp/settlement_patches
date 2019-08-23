PFPR=PfPR/Raster\ Data/PfPR_median/PfPR_median_Global_admin0_2017.tif
HRSL=hrsl/hrsl_uga_settlement.tif
HUMADMIN3=uga_admbnda/uga_admbnda_adm3_UBOS_v5.shp


all:  $(HRSL) $(PFPR) $(HUMADMIN3)

hrsl_uga_v1.zip:
	curl -o hrsl_uga_v1.zip https://www.ciesin.columbia.edu/repository/hrsl/hrsl_uga_v1.zip

$(HRSL): hrsl_uga_v1.zip
	mkdir -p hrsl
	(cd hrsl && unzip -n ../hrsl_uga_v1.zip)
	touch $(HRSL)

PfPR.zip:
	curl -o PfPR.zip https://map.ox.ac.uk/wp-content/uploads/lancet-2019/pf/PfPR.zip

$(PFPR): PfPR.zip
	-unzip PfPR.zip
	touch $(PFPR)

uga_admbnda_adm3_ubos_v5.zip:
	curl -o uga_admbnda_adm3_ubos_v5.zip https://data.humdata.org/dataset/6d6d1495-196b-49d0-86b9-dc9022cde8e7/resource/48de0e71-f1e8-4720-8e82-b60bfa36eab4/download/uga_admbnda_adm3_ubos_v5.zip
	curl -o uga_admbnda_adm2_ubos_v2.zip https://data.humdata.org/dataset/6d6d1495-196b-49d0-86b9-dc9022cde8e7/resource/5f36939d-cf46-499c-ba32-eec192239459/download/uga_admbnda_adm2_ubos_v2.zip
	curl -o uga_admbnda_adm1_ubos_v2.zip https://data.humdata.org/dataset/6d6d1495-196b-49d0-86b9-dc9022cde8e7/resource/207c13a1-7cd7-4045-9763-afed253b423a/download/uga_admbnda_adm1_ubos_v2.zip
	curl -o uga_admbnda_adm0_ubos_v2.zip https://data.humdata.org/dataset/6d6d1495-196b-49d0-86b9-dc9022cde8e7/resource/db0c97bb-c3cc-4b7c-908a-982cb5a6d881/download/uga_admbnda_adm0_ubos_v2.zip

$(HUMADMIN3): uga_admbnda_adm3_ubos_v5.zip
	mkdir -p uga_admbnda
	(cd uga_admbnda && ls ../uga_admbnda_adm*.zip | xargs -n 1 unzip)

uga_admin_energy.zip:
	curl -o uga_admin_energy.zip http://data.energy-gis.opendata.arcgis.com/datasets/f0d63758fb8f4ded85394b51594d294a_0.zip



clean:
	rm -rf hrsl PfPR *.zip
