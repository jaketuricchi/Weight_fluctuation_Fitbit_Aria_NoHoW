library('dplyr')
library('lubridate')
library('birk')
library('broom')
library('imputeMulti')
library('imputeTS')
library('zoo')
library('xts')
library('reshape2')
library('plyr')
library('Metrics')
library('Hmisc')
library('forecast')
library('tidyr')
library('psych')
library('forcats')
library('anytime')
library('rowr')
library('tableone')
library('tidyselect')
library('magrittr')
library('gridExtra')
library(ggplot2)
library('car')
library('multcomp')
options(scipen = 999)

#setwd
setwd("C:/Users/jaket/Dropbox/PhD/NoHoW Analyses/Weight variability")
setwd("~/Dropbox/PhD/NoHoW Analyses/Weight variability")

#read weight df in from aria smart scale data (already preprocessed with outliers removed)
df<-read.csv('Daily_data_with_NAs_020719.csv')%>%dplyr::select(ID, tseq, day_no, weight)%>% 
  filter(day_no<731)#2 years
colnames(df)[2]<-'date'

#dates to date
df$date<-anydate(df$date)

#find n_weights provided by each participant in order to set minimum criteria
df<-df%>%group_by(ID)%>%filter(!is.na(weight))%>%dplyr::summarise(n_weights=n())%>%
  merge(df, ., by='ID')

#C:/Users/jaket/iCloudDrive/
#~/Library/Mobile Documents/com~apple~CloudDocs/
#add participant characteristics (grouping variables, age, gender, centre, BMI) for description by group at each stage
characteristics  <-readRDS(file.choose())%>%
  dplyr::select(participant_study_id, elig_gender, elig_age, ecid1_height_recorded, ecid1_weight_recorded)%>%
  mutate(BMI=ecid1_weight_recorded/(ecid1_height_recorded/100)^2, 
         centre=as.factor(ifelse(grepl("UL_", participant_study_id, ignore.case = T), "UL", 
                                 ifelse(grepl("CPH_", participant_study_id, ignore.case = T), 
                                        "CPH", "LIS"))),
         age_group=as.factor(ifelse(elig_age<30, 'under 30', 
                                     ifelse(elig_age>=30 & elig_age <46, '30 to 45',
                                            ifelse(elig_age>=46 & elig_age <61, '45 to 60', 'over 60')))),
         BMI_group=as.factor(ifelse(BMI<25, 'Healthy_weight', 
                          ifelse(BMI>=25 & BMI <30, 'Overweight',
                                 ifelse(BMI>=30 & BMI <35, 'Obese_C1', 'Obese_C2_3')))))%>%
  dplyr::select(-ecid1_height_recorded)
colnames(characteristics)<-c('ID', 'Gender', 'Age', 'First_weight', 'BMI', 'Centre', 'AgeGroup', 'BMI_group')

##merge characteristics into weights df
df2<-merge(df, characteristics, by='ID')
df=NULL #remove these as we go as to not have several large dfs loaded.

#minimum days >= 20 for detrending
df3<-filter(df2, n_weights>19)
df2=NULL #remove this line if we need the unfiltered data

###############################
# day of week (dow) analysis #
###############################
dow_df<-df3

#to detrend to data we fit a loess. detrending accounts for weight change confounding fluctuation patterns
## the span of the loess will be optimised using a minimization fn
calcSSE <- function(x){
  loessMod <- try(loess(weight ~ day_no, data=ppt, span=x), silent=T)
  res <- try(loessMod$residuals, silent=T)
  if(class(res)!="try-error"){
    if((sum(res, na.rm=T) > 0)){
      sse <- sum(res^2)  
    }
  }else{
    sse <- 99999
  }
  return(sse)
}

#generate a fn to fit a loess to each ppts body weights and extract the relative residuals.
loess_fn<-function(x){
  span<-optim(par=c(0.5), calcSSE, method="CG")
  loess_weights=predict(loess(weight~day_no, data=x), x$day_no, span=span$par, degree=2)
  detrended_kg= x$weight - loess_weights
  x$relative_residual=(detrended_kg/x$weight)*100 #now we have the % residual relative to the polynomial trend
  return(x)
}

#apply loess fn
dow_df<-dlply(dow_df, 'ID', loess_fn)%>%bind_rows()

#add day of the week
dow_df$dow<-as.factor(weekdays(dow_df$date))

#define minimum criteria for DOW analysis - at least 1 weight reading on all 7 days at some given point
# (loose criteria for now to maximise retention)

eligible_dow<-dow_df%>%group_by(dow, ID)%>%dplyr::summarise(n_per_dow=n())%>%group_by(ID)%>%
  dplyr::summarise(days=n())%>%filter(days==7)

dow_df<-subset(dow_df, ID %in% eligible_dow$ID)

#sort factor levels
dow_df<-dow_df%>%mutate(dow=fct_relevel(dow,  "Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
                                        "Saturday", "Sunday"),
                        AgeGroup=fct_relevel(AgeGroup, "under 30", "30 to 45", "45 to 60", "over 60"),
                        BMI_group=fct_relevel(BMI_group, "Healthy_weight", "Overweight", "Obese_C1", "Obese_C2_3"))

#shorten day names for plots                        
dow_df$dow<-revalue(dow_df$dow, c("Monday"= "Mon", "Tuesday"="Tues", "Wednesday"="Wed", "Thursday"="Thurs", 
                                  "Friday"= "Fri", "Saturday"="Sat", "Sunday"= "Sun"))   

#calculate means for each day, keep seperate dfs as to not delete factor structures
dow_summary_all<-dow_df%>%group_by(dow)%>%dplyr::summarise(fluctuation=mean(relative_residual, na.rm=T),
                                                           sd=sd(relative_residual, na.rm=T), Gender='All')
dow_summary_gender<-dow_df%>%group_by(dow, Gender)%>%dplyr::summarise(fluctuation=mean(relative_residual, na.rm=T),
                                                                      sd=sd(relative_residual, na.rm=T))%>%
  bind_rows(., dow_summary_all) #we'll plot all and gender together
dow_summary_centre<-dow_df%>%group_by(dow, Centre)%>%dplyr::summarise(fluctuation=mean(relative_residual, na.rm=T),
                                                                      sd=sd(relative_residual, na.rm=T))
dow_summary_BMI<-dow_df%>%group_by(dow, BMI_group)%>%dplyr::summarise(fluctuation=mean(relative_residual, na.rm=T),
                                                                      sd=sd(relative_residual, na.rm=T))
dow_summary_age<-dow_df%>%group_by(dow, AgeGroup)%>%dplyr::summarise(fluctuation=mean(relative_residual, na.rm=T),
                                                                     sd=sd(relative_residual, na.rm=T))

#do plots
dow_gender_plot<-ggplot(dow_summary_gender, aes(x=dow, y=fluctuation, color=Gender, group=Gender))+
  geom_point(size=2.5)+geom_line(size=1)+
  ylab("Detrended weight (%)")+ xlab("Day of the week")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'Gender')

dow_centre_plot<-ggplot(dow_summary_centre, aes(x=dow, y=fluctuation, color=Centre, group=Centre))+
  geom_point(size=2.5)+geom_line(size=1)+
  ylab("Detrended weight (%)")+ xlab("Day of the week")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'Gender')

dow_BMI_plot<-ggplot(dow_summary_BMI, aes(x=dow, y=fluctuation, color=BMI_group, group=BMI_group))+
  geom_point(size=2.5)+geom_line(size=1)+
  ylab("Detrended weight (%)")+ xlab("Day of the week")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'Gender')

dow_age_plot<-ggplot(dow_summary_age, aes(x=dow, y=fluctuation, color=AgeGroup, group=AgeGroup))+
  geom_point(size=2.5)+geom_line(size=1)+
  ylab("Detrended weight (%)")+ xlab("Day of the week")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  scale_x_discrete(limits = c("Mon",'Tues','Wed','Thu', 'Fri', 'Sat', 'Sun'))+ 
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'Gender')

grid.arrange(dow_gender_plot, dow_centre_plot, dow_BMI_plot, dow_age_plot)

#####################
# group differences #
#####################
#here we test for group differences at each timepoint using Anova SS3.

#type 3 SS Anova function, flags significant values so we can post-hoc those
SS3_dow_fn<-function(x){
  SS3<-Anova(lm(relative_residual~Gender + Centre + BMI_group + AgeGroup, data=x), type='3')
  SS3<-data.frame(SS3)%>%mutate(ph=ifelse(Pr..F. < 0.05,1,0)) #ph = needs post-hoc test to probe group diffs.
  SS3$grouping<-rownames(SS3)
  SS3<-SS3[-1,]
  SS3$dow<-x$dow[1]
  return(SS3)
}

SS3_dow_results<-dlply(dow_df, 'dow', SS3_dow_fn)%>%bind_rows()
#write.csv(SS3_dow_results, 'SS3_dow_results.csv', row.names = F)


#function to extract tabulated Tukey post hoc results
get_tukey_results<-function(t){
  PH_df<-summary(t)$test
  PH_df2<-cbind(PH_df$coefficients, PH_df$sigma, PH_df$tstat,
                       PH_df$pvalues)
  error <- attr(PH_df$pvalues, "error")
  colnames(PH_df2) <- c("Estimate", "Std.Error","t.value", "p.value")
  PH_df2<-data.frame(round(PH_df2, 3))
  PH_df2$comparison<-rownames(PH_df2)
  return(PH_df2)
}

#run post-hocs extract results for group differences between each day
get_dow_anova_group_results<-function(x){
  aov_dow<-aov(relative_residual~Gender+ Centre + BMI_group + AgeGroup, data=x)
  PH_gender <- glht(aov_dow, linfct = mcp(Gender = "Tukey"))
  PH_centre <- glht(aov_dow, linfct = mcp(Centre = "Tukey"))
  PH_BMI_status<-glht(aov_dow, linfct = mcp(BMI_group = "Tukey"))
  PH_age_group<-glht(aov_dow, linfct = mcp(AgeGroup = "Tukey"))
  
  dow_tukey_gender_results<-get_tukey_results(PH_gender)%>%mutate(grouping='Gender')
  dow_tukey_centre_results<-get_tukey_results(PH_centre)%>%mutate(grouping='Centre')
  dow_tukey_BMI_results<-get_tukey_results(PH_BMI_status)%>%mutate(grouping='BMI_group')
  dow_tukey_age_results<-get_tukey_results(PH_age_group)%>%mutate(grouping='AgeGroup')
  
  dow_all_tukey_results<-bind_rows(dow_tukey_gender_results, dow_tukey_centre_results,
                                   dow_tukey_BMI_results, dow_tukey_age_results)%>%
    mutate(dow=x$dow[1])
  return(dow_all_tukey_results)
}

all_tukey_results_df<-dlply(dow_df, 'dow', get_dow_anova_group_results)%>%bind_rows() #this has all Tukey results

#We only want the Tukey results which were significant in the type3 Anova.
# so we need to merge and filter with the SS3 df

all_dow_anova_results<-merge(SS3_dow_results, all_tukey_results_df, by=c('dow', 'grouping'))

dow_tukey_significant<-all_dow_anova_results%>%filter(ph==1) #a list of all required post_hoc results
#write.csv(dow_tukey_significant, 'Post-hoc_results_dow.csv', row.names=F)

# DOW analyses complete.







######################
## holiday analysis ##
######################
# how much does body weight fluctuate upwards over the Christmas period after accounting for the wider trend?

#extract a wider period so we can get a trend over xmas, start by making all dates the same year to make things easier
holiday_df<-df3%>%
  mutate(inc=ifelse(date> '2017-11-27' & date < '2018-02-28', 1,
                        ifelse(date>'2018-11-27'& date < '2019-02-28', 2, NA)),
         year=as.numeric(format(date,'%Y')))%>%
  filter(!is.na(inc))%>%arrange(ID, day_no)

#we need to arbitrarily defin inclusion criteria. For this, we will consider this as
#2 weights before, during and after each christmas 'period' for inclusion.
#participants must meet this at at least 1 year.

elig18<-holiday_df%>%group_by(ID)%>%filter(!is.na(weight))%>%
  mutate(pre2018=ifelse(date>'2017-11-27' & date<'2017-12-18', 1, 0),
         during2018=ifelse(date>'2017-12-18' & date<'2018-01-08', 1, 0),
         post2018=ifelse(date>'2018-01-08' & date<'2018-02-28', 1, 0))%>%
  dplyr::summarise(total_pre=sum(pre2018), total_during=sum(during2018), total_post=sum(post2018))%>%
  filter(total_pre>= 2 & total_during>= 2 & total_post>= 2)
elig18$inc<-1

elig19<-holiday_df%>%group_by(ID)%>%filter(!is.na(weight))%>%
  mutate(pre2019=ifelse(date>'2018-11-27' & date<'2018-12-18', 1, 0),
         during2019=ifelse(date>'2018-12-18' & date<'2019-01-08', 1, 0),
         post2019=ifelse(date>'2019-01-08' & date<'2019-02-28', 1, 0))%>%
  dplyr::summarise(total_pre=sum(pre2019), total_during=sum(during2019), total_post=sum(post2019))%>%
  filter(total_pre>= 2, total_during>= 2, total_post>= 2)
elig19$inc<-2

eligible_holiday<-bind_rows(elig18, elig19)%>%dplyr::select(ID, inc) # these IDs are a list of people with enough data for holiday analysis

#impute weight using a wide EWMA since we are looking at wider fluctuation and not day-to-day variability
#split by id and year when fitting MA so that we fit a seperate MA for each year.
x<-filter(holiday_df, ID=='CPH_AKJ_8970')
ma_function_weight<-function(x){
  x<-x%>%group_by(inc)%>%
    mutate(weight_ma=na.ma(x$weight, k=3, weighting='exponential'))
  return(x)
}

holiday_df2<-merge(holiday_df, eligible_holiday, by=c('ID', 'inc'))%>%
  dlply(., 'ID', ma_function_weight)%>%bind_rows()
holiday_df=NULL

#now we need to fit a loess with a large span to the curve (large span = more linear)
holiday_loess_fn<-function(x){
    x$loess_weights=predict(loess(weight~day_no, data=x), x$day_no, span=0.8, degree=2)
    detrended_kg= x$weight_ma - x$loess_weights
    x$relative_residual=(detrended_kg/x$weight)*100 #now we have the % residual relative to the polynomial trend
    return(x)
  }

holiday_df3<-dlply(holiday_df2, 'ID', holiday_loess_fn)%>%bind_rows()
holiday_df2=NULL


#######################
#merge all years together, in order to keep date format this seems to be most suitable:
holiday_df3$year<-as.numeric(format(holiday_df3$date,'%Y'))
y2018<-holiday_df3%>%
  filter(year=='2018', date<'2018-04-01')%>%mutate(date_all=date)
y2018_2<-holiday_df3%>%
  filter(year=='2018', date>'2018-04-01')%>%mutate(date_all=date-years(1))
y2019<-holiday_df3%>%
  filter(year=='2019', date<'2019-04-01')%>%mutate(date_all=date-years(1))
y2017<-holiday_df3%>%
  filter(year=='2017')%>%mutate(date_all=date)

holiday_df4<-bind_rows(y2017, y2018, y2018_2, y2019)
holiday_df3=NULL

#summarise residuals by group for plotting
holiday_summary_all<-holiday_df4%>%group_by(date_all)%>%dplyr::summarise(fluctuation=mean(relative_residual, na.rm=T),
                                                           sd=sd(relative_residual, na.rm=T), Gender='All')
holiday_summary_gender<-holiday_df4%>%group_by(date_all, Gender)%>%dplyr::summarise(fluctuation=mean(relative_residual, na.rm=T),
                                                                      sd=sd(relative_residual, na.rm=T))%>%
  bind_rows(., holiday_summary_all) #we'll plot all and gender together
holiday_summary_centre<-holiday_df4%>%group_by(date_all, Centre)%>%dplyr::summarise(fluctuation=mean(relative_residual, na.rm=T),
                                                                      sd=sd(relative_residual, na.rm=T))
holiday_summary_BMI<-holiday_df4%>%group_by(date_all, BMI_group)%>%dplyr::summarise(fluctuation=mean(relative_residual, na.rm=T),
                                                                      sd=sd(relative_residual, na.rm=T))
holiday_summary_age<-holiday_df4%>%group_by(date_all, AgeGroup)%>%dplyr::summarise(fluctuation=mean(relative_residual, na.rm=T),
                                                                     sd=sd(relative_residual, na.rm=T))

#do plots
holiday_gender_plot<-ggplot(holiday_summary_gender, aes(x=date_all, y=fluctuation, color=Gender, group=Gender))+
  geom_point(size=2.5)+
  stat_smooth(method='loess', span=0.4, se=F)+
  ylab("Detrended weight (%)")+ xlab("Date")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'Gender')

holiday_centre_plot<-ggplot(holiday_summary_centre, aes(x=date_all, y=fluctuation, color=Centre, group=Centre))+
  geom_point(size=2.5)+
  stat_smooth(method='loess', span=0.4, se=F)+
  ylab("Detrended weight (%)")+ xlab("Date")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'Centre')

holiday_BMI_plot<-ggplot(holiday_summary_BMI, aes(x=date_all, y=fluctuation, color=BMI_group, group=BMI_group))+
  geom_point(size=2.5)+
  stat_smooth(method='loess', span=0.4, se=F)+
  ylab("Detrended weight (%)")+ xlab("Date")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'BMI_group')

holiday_age_plot<-ggplot(holiday_summary_age, aes(x=date_all, y=fluctuation, color=AgeGroup, group=AgeGroup))+
  geom_point(size=2.5)+
  stat_smooth(method='loess', span=0.4, se=F)+
  ylab("Detrended weight (%)")+ xlab("Date")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'AgeGroup')

grid.arrange(holiday_gender_plot, holiday_centre_plot, holiday_BMI_plot, holiday_age_plot)


#####################
# group differences #
#####################
# define holiday period gain in each person

holiday_gain_fn<-function(x){
  before<-x%>%group_by(inc)%>%filter(date_all < '2017-12-25')%>%
    dplyr::summarise(min_before=min(weight_ma))
  after<-x%>%group_by(inc)%>%filter(date_all > '2017-12-25')%>%
    dplyr::summarise(max_after=max(weight_ma))
  change<-merge(before, after)%>%mutate(difference=max_after - min_before,
                                        relative_diff=(difference/min_before)*100)
  change$ID<-x$ID[1]
  return(change)
}

#merge characteristics into weight gain so we can test group differences in gain.
holiday_weight_change<-dlply(holiday_df4, 'ID', holiday_gain_fn)%>%bind_rows()%>%
  merge(., characteristics, by='ID')%>%
  mutate(AgeGroup=fct_relevel(AgeGroup, "under 30", "30 to 45", "45 to 60", "over 60"),
         BMI_group=fct_relevel(BMI_group, "Healthy_weight", "Overweight", "Obese_C1", "Obese_C2_3"))

#mean and SD per group
holiday_gain_gender<-holiday_weight_change%>%group_by(Gender)%>%
  dplyr::summarise(mean_change=mean(relative_diff, na.rm=T), sd=sd(relative_diff, na.rm = T))
holiday_gain_centre<-holiday_weight_change%>%group_by(Centre)%>%
  dplyr::summarise(mean_change=mean(relative_diff, na.rm=T), sd=sd(relative_diff, na.rm = T))
holiday_gain_BMI<-holiday_weight_change%>%group_by(BMI_group)%>%
  dplyr::summarise(mean_change=mean(relative_diff, na.rm=T), sd=sd(relative_diff, na.rm = T))
holiday_gain_age<-holiday_weight_change%>%group_by(AgeGroup)%>%
  dplyr::summarise(mean_change=mean(relative_diff, na.rm=T), sd=sd(relative_diff, na.rm = T))

####### ANOVAs
#type 3 SS Anova , flags significant values so we can post-hoc those
SS3_results_holiday<-Anova(lm(relative_diff~Gender + Centre + BMI_group + AgeGroup, data=holiday_weight_change), type='3')
SS3_results_holiday<-data.frame(SS3_results_holiday)%>%
  mutate(ph=ifelse(Pr..F. < 0.05,1,0)) #ph = needs post-hoc test to probe group diffs.
SS3_results_holiday$grouping<-rownames(SS3_results_holiday)
#write.csv(SS3_results_holiday, 'Anova_results_holiday.csv', row.names=F)

### POSTHOC
aov_holiday<-aov(relative_diff~Gender+Centre+BMI_group+AgeGroup, data=holiday_weight_change)
PH_centre=glht(aov_holiday, linfct=mcp(Centre='Tukey'))
PH_centre_results<-get_tukey_results(PH_centre) # results for centre difference testing.
#write.csv(PH_centre_results, 'Post-hoc_results_holiday.csv', row.names=F)



#####################
# seasonal analysis #
#####################
# how does body weight fluctuate across seasons in different groups, after accounting for linear trends?

#first, define minimum criteria and filter
# 5 weights in each season are required.

seasonal_df<-df3%>%mutate(yearless_date=as.numeric(as.character(format(date, '%m%d'))),
                          season=ifelse(yearless_date >= 0320 & yearless_date < 0620, 'Spring',
                                           ifelse(yearless_date >= 0620 & yearless_date < 0922, 'Summer',
                                                  ifelse(yearless_date >= 0922 & yearless_date < 1222, 'Autumn',
                                                         ifelse(yearless_date >=1222 & yearless_date < 1231, 'Winter', 'Winter')))))

eligible_seasonal<-seasonal_df%>%group_by(ID, season)%>% filter(!is.na(weight))%>%
  dplyr::summarise(n=n())%>%group_by(ID)%>% filter(n>4)%>%dplyr::summarise(n=n())%>%
  filter(n==4) # this picks only IDs with at least 5 weights in all 4 seasons
  
seasonal_df2<-subset(seasonal_df, ID %in% eligible_seasonal$ID)
seasonal_df=NULL

#fit lm through the observed data
fit_lm_to_seasonal<-function(x){
  y<-na.omit(x)
  y$lm_weight<-predict(lm(weight~day_no, data=y))
  x2<-merge(x, y, all=T)
  x2$lm_weight<-na.interpolation(x2$lm_weight, option='linear')
  x2$weight_ma<-na.ma(x2$weight, k=3, weighting='exponential')
  x2$fluctuation_kg<-x2$weight_ma-x2$lm_weight
  x2$perc_fluctuation<-(x2$fluctuation_kg/x2$weight_ma)*100
  return(x2)
}

seasonal_df3<-dlply(seasonal_df2, 'ID', fit_lm_to_seasonal)%>%bind_rows()
seasonal_df2=NULL

#merge all onto 1 year
seasonal_df3$year<-as.numeric(format(seasonal_df3$date,'%Y'))
y2017<-seasonal_df3%>%
  filter(year=='2017')%>%mutate(date_all=date+years(1))
y2018<-seasonal_df3%>%
  filter(year=='2018')%>%mutate(date_all=date)
y2019<-seasonal_df3%>%
  filter(year=='2019')%>%mutate(date_all=date-years(1))
 
seasonal_df4<-bind_rows(y2017, y2018, y2019)
seasonal_df3=NULL 

#summarise data by day of year


#################
## plot seasons #
#################

seasonal_summary_all<-seasonal_df4%>%group_by(date_all)%>%dplyr::summarise(fluctuation=mean(perc_fluctuation, na.rm=T),
                                                                         sd=sd(perc_fluctuation, na.rm=T), Gender='All')
seasonal_summary_gender<-seasonal_df4%>%group_by(date_all, Gender)%>%dplyr::summarise(fluctuation=mean(perc_fluctuation, na.rm=T),
                                                                                    sd=sd(perc_fluctuation, na.rm=T))%>%
  bind_rows(., seasonal_summary_all) #we'll plot all and gender together
seasonal_summary_centre<-seasonal_df4%>%group_by(date_all, Centre)%>%dplyr::summarise(fluctuation=mean(perc_fluctuation, na.rm=T),
                                                                                    sd=sd(perc_fluctuation, na.rm=T))
seasonal_summary_BMI<-seasonal_df4%>%group_by(date_all, BMI_group)%>%dplyr::summarise(fluctuation=mean(perc_fluctuation, na.rm=T),
                                                                                    sd=sd(perc_fluctuation, na.rm=T))
seasonal_summary_age<-seasonal_df4%>%group_by(date_all, AgeGroup)%>%dplyr::summarise(fluctuation=mean(perc_fluctuation, na.rm=T),
                                                                                   sd=sd(perc_fluctuation, na.rm=T))

seasonal_gender_plot<-ggplot(seasonal_summary_gender, aes(x=date_all, y=fluctuation, color=Gender, group=Gender))+
  stat_smooth(method='loess', span=0.4, se=F)+
  ylab("Detrended weight (%)")+ xlab("Date")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'Gender')+
  geom_vline(xintercept = as.numeric(as.Date('20-03-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('20-06-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('22-09-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('21-12-2018', format='%d-%m-%Y')), linetype=2, size=1)

seasonal_centre_plot<-ggplot(seasonal_summary_centre, aes(x=date_all, y=fluctuation, color=Centre, group=Centre))+
  stat_smooth(method='loess', span=0.4, se=F)+
  ylab("Detrended weight (%)")+ xlab("Date")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'Centre')+
  geom_vline(xintercept = as.numeric(as.Date('20-03-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('20-06-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('22-09-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('21-12-2018', format='%d-%m-%Y')), linetype=2, size=1)

seasonal_BMI_plot<-ggplot(seasonal_summary_BMI, aes(x=date_all, y=fluctuation, color=BMI_group, group=BMI_group))+
  stat_smooth(method='loess', span=0.4, se=F)+
  ylab("Detrended weight (%)")+ xlab("Date")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'BMI_group')+
  geom_vline(xintercept = as.numeric(as.Date('20-03-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('20-06-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('22-09-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('21-12-2018', format='%d-%m-%Y')), linetype=2, size=1)

seasonal_age_plot<-ggplot(seasonal_summary_age, aes(x=date_all, y=fluctuation, color=AgeGroup, group=AgeGroup))+
  stat_smooth(method='loess', span=0.4, se=F)+
  ylab("Detrended weight (%)")+ xlab("Date")+
  theme(axis.title = element_text(size=18, color='black'))+ 
  theme(axis.text=element_text(size=16, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 18),legend.text = element_text(size=18))+
  scale_fill_discrete(name = 'AgeGroup')+
  geom_vline(xintercept = as.numeric(as.Date('20-03-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('20-06-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('22-09-2018', format='%d-%m-%Y')), linetype=2, size=1)+
  geom_vline(xintercept = as.numeric(as.Date('21-12-2018', format='%d-%m-%Y')), linetype=2, size=1)

grid.arrange(seasonal_gender_plot, seasonal_centre_plot, seasonal_BMI_plot, seasonal_age_plot)


####################
# group differences #
#####################

#get means for each season
seasonal_means<-seasonal_df4%>%group_by(ID, season)%>%dplyr::summarise(mean_season_fluctuation=mean(perc_fluctuation))%>%
  merge(., characteristics, by='ID')%>%
  mutate(AgeGroup=fct_relevel(AgeGroup, "under 30", "30 to 45", "45 to 60", "over 60"),
         BMI_group=fct_relevel(BMI_group, "Healthy_weight", "Overweight", "Obese_C1", "Obese_C2_3"))

#type 3 SS Anova function, flags significant values so we can post-hoc those
SS3_season_fn<-function(x){
  SS3<-Anova(lm(mean_season_fluctuation~Gender + Centre + BMI_group + AgeGroup, data=x), type='3')
  SS3<-data.frame(SS3)%>%mutate(ph=ifelse(Pr..F. < 0.05,1,0)) #ph = needs post-hoc test to probe group diffs.
  SS3$grouping<-rownames(SS3)
  SS3<-SS3[-1,]
  SS3$season<-x$season[1]
  return(SS3)
}

SS3_season_results<-dlply(seasonal_means, 'season', SS3_season_fn)%>%bind_rows()
#write.csv(SS3_season_results, 'SS3_season_results.csv', row.names = F)


#run post-hocs extract results for group differences between each day
get_season_anova_group_results<-function(x){
  aov_season<-aov(mean_season_fluctuation~Gender+ Centre + BMI_group + AgeGroup, data=x)
  PH_gender <- glht(aov_season, linfct = mcp(Gender = "Tukey"))
  PH_centre <- glht(aov_season, linfct = mcp(Centre = "Tukey"))
  PH_BMI_status<-glht(aov_season, linfct = mcp(BMI_group = "Tukey"))
  PH_age_group<-glht(aov_season, linfct = mcp(AgeGroup = "Tukey"))
  
  season_tukey_gender_results<-get_tukey_results(PH_gender)%>%mutate(grouping='Gender')
  season_tukey_centre_results<-get_tukey_results(PH_centre)%>%mutate(grouping='Centre')
  season_tukey_BMI_results<-get_tukey_results(PH_BMI_status)%>%mutate(grouping='BMI_group')
  season_tukey_age_results<-get_tukey_results(PH_age_group)%>%mutate(grouping='AgeGroup')
  
  season_all_tukey_results<-bind_rows(season_tukey_gender_results, season_tukey_centre_results,
                                   season_tukey_BMI_results, season_tukey_age_results)%>%
    mutate(season=x$season[1])
  return(season_all_tukey_results)
}

all_tukey_results_df<-dlply(seasonal_means, 'season', get_season_anova_group_results)%>%bind_rows() #this has all Tukey results

#We only want the Tukey results which were significant in the type3 Anova.
# so we need to merge and filter with the SS3 df

all_season_anova_results<-merge(SS3_season_results, all_tukey_results_df, by=c('season', 'grouping'))

season_tukey_significant<-all_season_anova_results%>%filter(ph==1) #a list of all required post_hoc results
write.csv(season_tukey_significant, 'Tukey_seasonal_results.csv', row.names = F)


######################################
# figures and tables for publication #
######################################

#sample characteristics - since there are 3 samples we will report this 3x.

#what was the mean measurement duration and mean scale use over this duration?
duration<-df3%>%group_by(ID)%>%dplyr::summarise(duration_days=n())
weights<-df3%>%group_by(ID)%>%filter(!is.na(weight))%>%dplyr::summarise(n_weights=n())

#merge these in with characteristics
descriptives<-merge(characteristics, duration)
descriptives2<-merge(descriptives, weights)

#first, get these 3 eligible samples together in 1 df and define the 3 samples using binary vars
dow_sample<-subset(df3, ID %in% dow_df$ID)
dow_sample$dow_sample=1
dow_sample$group='dow'

holiday_sample<-subset(df3, ID %in% holiday_df4$ID)
holiday_sample$holiday_sample=1
holiday_sample$group='holiday'

season_sample<-subset(df3, ID %in% seasonal_df4$ID)
season_sample$season_sample=1
season_sample$group='season'

all_samples<-merge(dow_sample, holiday_sample, all=T)
all_samples2<-merge(all_samples, season_sample, all=T)

descriptives3<-merge(descriptives2, all_samples2, all=T)%>%
  group_by(ID)%>%
  slice(1:1)

#get rid of big dfs
dow_sample=NULL
holiday_sample=NULL
season_sample=NULL
all_samples=NULL

descriptives3$First_weight<-as.numeric(descriptives3$First_weight)

##generate descriptive tables using TableOne
dow_descriptives<-CreateTableOne(vars = c('Gender', 'Centre', 'AgeGroup', 'BMI_Group',
                                          'First_weight', 'duration_days', 'n_weights'),
                                 data=filter(descriptives3, dow_sample==1))
print(dow_descriptives,showAllLevels = TRUE)
#write.csv(print(dow_descriptives, showAllLevels = TRUE), 
          #file='dow_sample_descriptives')

holiday_descriptives<-CreateTableOne(vars = c('Gender', 'Centre', 'AgeGroup', 'BMI_Group',
                                          'First_weight', 'duration_days', 'n_weights'),
                                 data=filter(descriptives3, holiday_sample==1))
print(holiday_descriptives,showAllLevels = TRUE)
#write.csv(print(holiday_descriptives, showAllLevels = TRUE), 
  #file='holiday_sample_descriptives')

season_descriptives<-CreateTableOne(vars = c('Gender', 'Centre', 'AgeGroup', 'BMI_Group',
                                              'First_weight', 'duration_days', 'n_weights'),
                                     data=filter(descriptives3, season_sample==1))
print(season_descriptives,showAllLevels = TRUE)
#write.csv(print(season_descriptives, showAllLevels = TRUE), 
#file='season_sample_descriptives')

#######################
# descriptive figures #
#######################
# Describing the frequency of scale use by day of the week, month of year, and
# with regards to days in trial.

all_samples3<-bind_rows(dow_sample, holiday_sample, season_sample)%>%
  mutate(month=as.factor(format(date, "%B")), day=weekdays(date), 
         day=fct_relevel(day,  "Monday", "Tuesday", "Wednesday", "Thursday", "Friday",
                             "Saturday", "Sunday"))

#determine the completeness (%) of data collection by month
month_plot_possible<-all_samples3%>%group_by(group, month)%>%
  dplyr::summarise(available=n()) #total number of days available
month_plot_available<-all_samples3%>%group_by(group, month)%>%
  filter(!is.na(weight))%>% dplyr::summarise(possible=n()) #number of weights present
month_plot<-merge(month_plot_possible, month_plot_available)%>%
  mutate(perc_complete_date=(possible/available)*100)

#plot barchat for each group
monthly_freq<-ggplot(month_plot, aes(x=month, y=perc_complete_date, fill=group, color=group))+
  geom_bar(stat="identity", position='dodge')+
  xlab("Month of the year")+ylab("Available data (%)")+ 
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18, color=))+ theme(axis.text=element_text(size=18))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text= element_text(color='black'))+
  coord_flip()+ scale_fill_discrete(name = 'Sample', labels = c('Daily', 'Seasonal', 'Holiday'))

###################################
# same again with day of the week #


#determine the completeness (%) of data collection by day
day_plot_possible<-all_samples3%>%group_by(group, day)%>%
  dplyr::summarise(available=n()) #total number of days available
day_plot_available<-all_samples3%>%group_by(group, day)%>%
  filter(!is.na(weight))%>% dplyr::summarise(possible=n()) #number of weights present
day_plot<-merge(day_plot_possible, day_plot_available)%>%
  mutate(perc_complete_date=(possible/available)*100)

#plot barchat for each group
daily_freq<-ggplot(day_plot, aes(x=day, y=perc_complete_date, fill=group, color=group))+
  geom_bar(stat="identity", position='dodge')+
  xlab("day of the year")+ylab("Available data (%)")+ 
  theme(axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18, color=))+ theme(axis.text=element_text(size=18))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text= element_text(color='black'))+
  coord_flip()+ scale_fill_discrete(name = 'Sample', labels = c('Daily', 'Seasonal', 'Holiday'))

######################
# scale use over time #
#######################


LT_use<-all_samples2%>%mutate(week_no=floor(day_no/7)+1)%>%
  group_by(ID, week_no)%>%
  dplyr::summarise(total=n(), available=sum(!is.na(weight)))%>%
  group_by(week_no)%>%
  dplyr::summarise(mean_week=mean(available), se_week=birk::se(available))

#LT_use
weekly_plot<-ggplot(subset(LT_use, week_no<74), aes(x=week_no, y=mean_week))+
  geom_point(size=2)+geom_line(size=1.5)+xlim(0,78)+
  theme(axis.title = element_text(size=20, color='black'))+ 
  theme(axis.text=element_text(size=19, color='black'))+ 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(legend.title = element_text(size = 20),legend.text = element_text(size=20))+
  xlab('Week of Trial')+ylab('Average scale use per week')+ylim(2,5)+
  geom_errorbar(aes(ymin=mean_week-se_week, ymax=mean_week+se_week), width=0.5)
weekly_plot

lay <- rbind(c(1,1),
             c(2,3))
grid.arrange(weekly_plot,daily_freq,daily_freq, layout_matrix = lay)














