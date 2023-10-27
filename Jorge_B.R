library(tidyr)
library(dplyr)
library(readxl)
library(ggplot2)
library(stringr)
library(forcats)
library(envalysis)
library(soiltexture)
library(plotrix)
library(broom)
library(ggthemes)
library(effectsize)


path<- file.path("tables_1.xlsx")
# tables_1 contains the measured variables from all the experiments
#The first sheet contains experiments Active carbon, penetration test, pH, and texture
#The second sheet contains experiments cec, and erosion resistance
#The third sheet contains experiments bulk density and soil moisture

#Nomenclature sample:
#       SA__NM___2___R2
#.      P1__P2___P3__P4 
# P1 Garden - SA = StadtAcker  ES = EssbareStadt  FW = Feuerwache  OI = O'pflanz ist   PP = Petuelpark
# P2: Management - MM = Managed site  NM = Non-Managed site
# P3: Depth category - 1 = 0-10 cm   2 = 10-20 cm
# P4: Repetition number - R1 = Repetition 1    R2 = Repetition 2 



var_1<- read_excel(path,sheet = 1) # Importing first sheet of tables_1
var_2<- read_excel(path,sheet = 2) # Importing second sheet of tables_1
var_3<- read_excel(path,sheet = 3) # Importing third sheet of tables_1

theme_set(theme_classic())


#Renaming variable full_name to have "_"
var_1<-var_1 %>% 
  mutate(full_name = str_replace_all(full_name, " ", "_"))
var_2<-var_2 %>%
  mutate(full_name = str_replace_all(full_name, " ","_")) 



#Transforming name elements into factors for better analysis
var_1<-var_1%>%
  mutate(garden = as_factor(garden),management = factor(management, ordered = TRUE, 
                                                        levels = c("MM","NM")),depth_cat = as_factor(depth_cat), rep = as_factor(rep), pen=as.numeric(pen))
var_2<-var_2%>%
  mutate(garden = as_factor(garden),management = factor(management, ordered = TRUE, 
                                                        levels = c("MM","NM")),depth_cat = as_factor(depth_cat), rep = as_factor(rep))
var_3<-var_3%>%
  mutate(garden = as_factor(garden),management = factor(management, ordered = TRUE, 
                                                        levels = c("MM","NM")),depth_cat = as_factor(depth_cat), rep = as_factor(rep),
         lid_weight =as.numeric(lid_weight), wet_weight = as.numeric(wet_weight), dry_weight = as.numeric(dry_weight))


#Create tibble with names
names<- var_1 %>% select(n,full_name)
names2<- var_1 %>% select(n,full_name, garden, management,depth_cat,rep)

#Adding sample number to var_3 for better processing
var_3<-var_3%>%
 inner_join(names, by = "full_name")%>%
  select(n, full_name, garden, management, depth_cat, rep,subsamp,cyl_num,cyl_weight,lid_weight,wet_weight,dry_weight)


#Calculation of Bulk density and Moisture content as bd_m
bd_m<-var_3%>%
 mutate(cyl_comp = cyl_weight - lid_weight, cyl_vol = 100, dry_soil= dry_weight - cyl_comp, wet_soil = wet_weight - cyl_comp, 
        bulk_density = dry_soil/cyl_vol)%>%
  mutate (moisture_content = 100 * ((wet_soil-dry_soil)/dry_soil))%>%
  group_by(n)%>%
  summarize(mean_bd = mean(bulk_density), mean_moist = mean (moisture_content))%>%
  inner_join(names2, by="n")%>%
  select(n, mean_bd,mean_moist, management, garden,depth_cat,rep)

#Calculation of Active Carbon
ac<-var_1%>%
  select(n, garden, management,depth_cat, rep,weight_ac, abs_ac)%>%
  mutate(mol_final= 0.0499*abs_ac + 0.0005)%>%#Using the slope and intercept obtained from the calibration 
  mutate(mol_diff = 0.02 - mol_final, equiv = 9000, soil_conc = 0.02*1000/weight_ac, active_carbon = mol_diff * equiv * soil_conc)

#Calculation of CEC
cec<-var_2%>%
  select(n, garden, management,depth_cat, rep,cec_weight,abs_cec)%>%
  mutate(c_1 = 44.321*abs_cec - 0.6262, c_2 = c_1*4, n_1 = c_2*0.2,  cec = n_1 / 10*cec_weight) 


texture<-var_1%>%
  select(n, garden, management, depth_cat, rep, texture_weight, sand_weight, silt_weight, clay_weight)%>%
  mutate(SAND= 100*sand_weight/texture_weight, SILT = 100*silt_weight/texture_weight, CLAY = 100*clay_weight/texture_weight)


#Creation of data frame for plotting soil textures of managed sites
texture1<-texture%>%
  filter(management == "MM")%>%
  select(SAND, SILT, CLAY)

#Creation of data frame for plotting soil textures of unmanaged sites
texture2<-texture%>%
  filter(management == "NM")%>%
  select(SAND, SILT, CLAY)

#Creation of data frame for classification of soil textures of unmanaged sites
texture3<-var_1%>%
  select(n, garden, management, depth_cat, rep, texture_weight, sand_weight, silt_weight, clay_weight)%>%
  mutate(SAND= 100*sand_weight/texture_weight, SILT = 100*silt_weight/texture_weight, CLAY = 100*clay_weight/texture_weight)%>%
  filter(management == "NM")%>%
  select(CLAY, SILT, SAND)%>%#important to keep order of clay silt and sand for the soil texture package
  as.data.frame() #important for the function that the data is a data frame and not a tibble

#Creation of data frame for classification of soil textures of managed sites
texture4<-var_1%>%
  select(n, garden, management, depth_cat, rep, texture_weight, sand_weight, silt_weight, clay_weight)%>%
  mutate(SAND= 100*sand_weight/texture_weight, SILT = 100*silt_weight/texture_weight, CLAY = 100*clay_weight/texture_weight)%>%
  filter(management == "MM")%>%
  select(CLAY, SILT, SAND)%>%#important to keep order of clay silt and sand for the soil texture package
  as.data.frame() #important for the function that the data is a data frame and not a tibble


#Importing soil classification

soil_cat <- as_tibble(TT.classes.tbl( class.sys = "USDA.TT" ))%>%
  select(abbr,name)%>%
  mutate(abbr = factor(abbr, levels = abbr))
soil_cat_fac<- soil_cat$abbr

#Obtaining soil texture classification  for unmanaged sites

res_textnm<-as_tibble(TT.points.in.classes(
  tri.data =texture3,
  class.sys = "USDA.TT"))%>%
  mutate(Cl = case_when( Cl == 1 ~ "Cl", Cl == 0 ~ NA),
         SiCl = case_when( SiCl == 1 ~ "SiCl", SiCl == 0 ~ NA),
         SaCl = case_when( SaCl == 1 ~ "SaCl", SaCl == 0 ~ NA),
         ClLo = case_when( ClLo == 1 ~ "ClLo", ClLo == 0 ~ NA),
         SiClLo  = case_when( SiClLo  == 1 ~ "SiClLo", SiClLo  == 0 ~ NA),
         SaClLo  = case_when( SaClLo  == 1 ~ "SaClLo", SaClLo  == 0 ~ NA),
         Lo  = case_when( Lo  == 1 ~ "Lo", Lo  == 0 ~ NA),
         SiLo  = case_when( SiLo  == 1 ~ "SiLo", SiLo  == 0 ~ NA),
         SaLo  = case_when( SaLo  == 1 ~ "SaLo", SaLo  == 0 ~ NA),
         Si  = case_when( Si  == 1 ~ "Si", Si  == 0 ~ NA),
         LoSa  = case_when( LoSa  == 1 ~ "LoSa", LoSa  == 0 ~ NA),
         Sa  = case_when( Sa  == 1 ~ "Sa", Sa  == 0 ~ NA))%>%
  unite("texture_cat",Cl:Sa, remove=TRUE, na.rm=TRUE)

#Obtaining soil texture classification  for managed sites

res_textmm<-as_tibble(TT.points.in.classes(
  tri.data =texture4,
  class.sys = "USDA.TT"))%>%
  mutate(Cl = case_when( Cl == 1 ~ "Cl", Cl == 0 ~ NA),
         SiCl = case_when( SiCl == 1 ~ "SiCl", SiCl == 0 ~ NA),
         SaCl = case_when( SaCl == 1 ~ "SaCl", SaCl == 0 ~ NA),
         ClLo = case_when( ClLo == 1 ~ "ClLo", ClLo == 0 ~ NA),
         SiClLo  = case_when( SiClLo  == 1 ~ "SiClLo", SiClLo  == 0 ~ NA),
         SaClLo  = case_when( SaClLo  == 1 ~ "SaClLo", SaClLo  == 0 ~ NA),
         Lo  = case_when( Lo  == 1 ~ "Lo", Lo  == 0 ~ NA),
         SiLo  = case_when( SiLo  == 1 ~ "SiLo", SiLo  == 0 ~ NA),
         SaLo  = case_when( SaLo  == 1 ~ "SaLo", SaLo  == 0 ~ NA),
         Si  = case_when( Si  == 1 ~ "Si", Si  == 0 ~ NA),
         LoSa  = case_when( LoSa  == 1 ~ "LoSa", LoSa  == 0 ~ NA),
         Sa  = case_when( Sa  == 1 ~ "Sa", Sa  == 0 ~ NA))%>%
  unite("texture_cat",Cl:Sa, remove=TRUE, na.rm=TRUE)


#Final matrix with texture categories.
  
texture_finnm<-var_1%>%
  select(n, garden, management, depth_cat, rep, texture_weight, sand_weight, silt_weight, clay_weight)%>%
  filter(management == "NM")%>%
  mutate(texture_class = factor(res_textnm$texture_cat,levels = soil_cat_fac))%>%
  left_join(soil_cat, by = c("texture_class" = "abbr"))%>%
  rename(text_class_long = name)

texture_finmm<-var_1%>%
  select(n, garden, management, depth_cat, rep, texture_weight, sand_weight, silt_weight, clay_weight)%>%
  filter(management == "MM")%>%
  mutate(texture_class = factor(res_textmm$texture_cat,levels = soil_cat_fac))%>%
  left_join(soil_cat, by = c("texture_class" = "abbr"))%>%
  rename(text_class_long = name)

#Calculation of water stable aggregates

water_stable<-var_2%>%
  select(n, full_name, garden, management, depth_cat, rep, soil_er_weight, g250_weight, recovery)%>%
  mutate(sand_frac = texture$SAND/100, sand_sample =  g250_weight*sand_frac, stable_agg = g250_weight - sand_sample, WSA = stable_agg/(recovery-sand_sample))%>%
    select(-sand_frac)

depth_label<- c("Depth: 0-10 cm", "Depth: 10-20 cm")
names(depth_label) <- c("1", "2")



# Box plots for all samples
 
 plot1b<-bd_m%>%
   ggplot(aes(management,mean_bd))+
   geom_boxplot()+
   labs(x = "Management", y = bquote('Bulk density '(g/cm^3)) )+
   scale_x_discrete(labels = c('Managed','Non-managed'))
 
 plot2b<-ac%>%
   ggplot(aes(management, active_carbon))+
   geom_boxplot()+
   labs(x = "Management", y = "Active carbon (mg/kg)")+
   scale_x_discrete(labels = c('Managed','Non-managed'))
 
 plot3b<-var_1%>%
   ggplot(aes(management,ph))+
   geom_boxplot()+
   labs(x = "Management", y = "pH")+
   scale_x_discrete(labels = c('Managed','Non-managed'))
 
 plot4b<- var_1%>%
   ggplot(aes(management,pen))+
   geom_boxplot()+
   labs(x = "Management", y = bquote("Soil hardness "(kg/cm^2)))+
   scale_x_discrete(labels = c('Managed','Non-managed'))
 
 
 plot5b<- cec%>%
   ggplot(aes(management,cec))+
   geom_boxplot()+
   labs(x = "Management", y = "Cation exchange capacity (cmol/kg)")+
   scale_x_discrete(labels = c('Managed','Non-managed'))
 
 plot6b<-bd_m%>%
   ggplot(aes(management,mean_moist))+
   geom_boxplot()+
   labs(x = "Management", y = "Moisture content (%)")+
   scale_x_discrete(labels = c('Managed','Non-managed'))
 
 
 plot8b<-water_stable%>%
   ggplot(aes(management, WSA))+
   geom_boxplot()+
   labs(x = "Management", y = bquote("Water stable aggregates " (g/g[soil])))+
   scale_x_discrete(labels = c('Managed','Non-managed'))
 
 
 
# Creation of variable with column names
 
test_name<-c("Bulk density","Soil moisture","Active carbon","CEC","Aggregate stability","Soil hardness","pH")


#t test results for all variables

t_results_all<-tidy(t.test(mean_bd ~ management, data = bd_m))%>%
  rbind(tidy(t.test(mean_moist ~ management, data = bd_m)))%>%
  rbind(tidy(t.test(active_carbon ~ management, data = ac)))%>%
  rbind(tidy(t.test(cec ~ management, data = cec)))%>%
  rbind(tidy(t.test(WSA ~ management, data = water_stable)))%>%
  rbind(tidy(t.test(pen ~ management, data = var_1)))%>%
  rbind(tidy(t.test(ph ~ management, data = var_1)))%>%
  transmute(test = test_name, difference = estimate, median_MM = estimate1, median_NM = estimate2, t= statistic, p_value = p.value )



#MANOVA test for assessing significance of difference between Managed and non-managed sites

full_data <- bd_m %>% select (mean_bd, mean_moist) %>% 
  cbind(ac$active_carbon)%>%
  cbind(cec$cec)%>%
  cbind(water_stable$WSA)%>%
  cbind(var_1$ph)%>%
  cbind(var_1$pen)%>%
  cbind(ac$management)%>%
  as_tibble()

names(full_data)<-c("bd","moist","ac","cec","WSA","ph","pen","management")
full_data<-full_data %>% mutate(management = as.character(management))
full_data[is.na(full_data)]<-0

dependent<-cbind(full_data$bd,full_data$moist,full_data$ac,full_data$cec,full_data$WSA,full_data$ph,full_data$pen)
manova_model<-manova(dependent ~ management,data = full_data)
summary(manova_model)
eta_squared(manova_model)


dependent_phys<-cbind(full_data$bd,full_data$moist,full_data$WSA,full_data$pen)
manova_model1<-manova(dependent_phys ~ management,data = full_data)
summary(manova_model1)
eta_squared(manova_model1)

dependent_chem<-cbind(full_data$ph,full_data$cec,full_data$ac)
manova_model2<-manova(dependent_chem ~ management,data = full_data)
summary(manova_model2)
eta_squared(manova_model2)
