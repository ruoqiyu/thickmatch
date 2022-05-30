\name{nysr}
\alias{nysr}
\docType{data}
\title{
Adolescent Work Intensity and Substance Use
}
\description{
NYSR data on adolescent work intensity and substance Use.
}
\usage{data("nysr")}
\format{
  A data frame with 2816 observations on the following 18 variables.
  \describe{
    \item{\code{IDS}}{NYSR identification number}
    \item{\code{intense}}{Based on question ``During the school year, about how many hours per week did you normally work at a paid job, or did you not have a job".
``Never": student did not have a job;
``Moderate": 1-19 hours;
``Intense": >=20 hours.
}
    \item{\code{family.income}}{Household income with 5000 = (between 0-10,000), 15000= (between 10,000 and 20,000),..., 95000 = (between 90,000 and 100,000) and 105,000 (above 100,000).}
    \item{\code{family.income.impute}}{Household income with 5000 = (between 0-10,000), 15000= (between 10,000 and 20,000),..., 95000 = (between 90,000 and 100,000) and 105,000 (above 100,000).  For subjects with missing family income, the mean is imputed.}
    \item{\code{family.income.mis}}{dummy variable for whether household income is missing and the mean is imputed.}
    \item{\code{family.structure}}{``Two Parent Biological": both biological father and mother living with child; ``Two Parent Nonbiological": someone assuming a mother role (biological, adoptive, stepparent) living with a husband who assumes a father role (biological, adoptive, step parent) where both parents are biological; ``Single Parent/Other": any other living situation for child.
}
    \item{\code{highest.education.parent.in.household}}{Maximum education level of household resident who assumes a mother role (biological, adoptive, stepparent) and household resident who assumes a father role (biological, adoptive, stepparent).  If the child is living with a single parent, then this is just the education level of that single parent.}
    \item{\code{highest.education.parent.in.household.impute}}{Maximum education level of household resident who assumes a mother role (biological, adoptive, stepparent) and household resident who assumes a father role (biological, adoptive, stepparent).  If the child is living with a single parent, then this is just the education level of that single parent. For subjects with missing highest education of parent in household, the mean is imputed.}
    \item{\code{highest.education.parent.in.household.mis}}{Dummy variable for whether household income is missing and the mean is imputed.}
    \item{\code{female}}{1 = female, 0 = male}
    \item{\code{race.black}}{1=black race, 0=other}
    \item{\code{race.hispanic}}{1=hispanic race, 0=other}
    \item{\code{age.teenager}}{age of teenager.  Age is imputed with the mean if it is missing}
    \item{\code{school.dropout}}{Dummy variable of whether student has dropped out of school}
    \item{\code{alcohol.use}}{Based on question ``How often, if at all, do you drink alcohol, such as beer, wine or mixed drinks, not including at religious services".
``Never": answered ``Never";
``Moderate": answered ``A few times year" or ``About once a month";
``Heavy": answered ``A few times a month", ``About once a week", ``A few times a week" or ``Almost every day".
}
    \item{\code{marijuana.use}}{Based on question ``How often, if ever, have you used marijuana?".
``Never": answered ``Never";
``Experimenter" answered ``You tried it once or twice";
``Continual User": answered ``You use it occasionally" or ``You use it regularly".
}
    \item{\code{p}}{Propensity score.
}
    \item{\code{plogit}}{Logit of propensity score.
}
  }
}
\details{
The following code constructed the data as used here.
wave1data$family.income=rep(NA,nrow(wave1data))
wave1data$family.income[wave1data$PINCOME1==1 & wave1data$PINCOME2==4]=5000
wave1data$family.income[wave1data$PINCOME1==1 & wave1data$PINCOME2==3]=15000
wave1data$family.income[wave1data$PINCOME1==1 & wave1data$PINCOME2==2]=25000
wave1data$family.income[wave1data$PINCOME1==1 & wave1data$PINCOME2==1]=35000
wave1data$family.income[wave1data$PINCOME1==2 & wave1data$PINCOME3==1]=45000
wave1data$family.income[wave1data$PINCOME1==2 & wave1data$PINCOME3==2]=55000
wave1data$family.income[wave1data$PINCOME1==2 & wave1data$PINCOME3==3]=65000
wave1data$family.income[wave1data$PINCOME1==2 & wave1data$PINCOME3==4]=75000
wave1data$family.income[wave1data$PINCOME1==2 & wave1data$PINCOME3==5]=85000
wave1data$family.income[wave1data$PINCOME1==2 & wave1data$PINCOME3==6]=95000
wave1data$family.income[wave1data$PINCOME1==2 & wave1data$PINCOME3==7]=105000
# For subjects with missing family income data, fill in mean and create a missing data indicator
wave1data$family.income.mis=is.na(wave1data$family.income)
#wave1data$family.income[wave1data$family.income.mis==1]=mean(wave1data$family.income,na.rm=TRUE)

# Find family structure variable
wave1data$family.structure=rep(NA,nrow(wave1data))
# wave1data$family.structure[wave1data$PMOTHER==1 & wave1data$PLIVE==1 & wave1data$PSPRELAT==1]="Two Parent Biological"
# wave1data$family.structure[wave1data$PMOTHER==1 & wave1data$PLIVE==2 & wave1data$PPARTPAR==1]="Two Parent Biological"
# wave1data$family.structure[wave1data$PFATHER==1 & wave1data$PLIVE==1 & wave1data$PSPRELAT==1]="Two Parent Biological"
# wave1data$family.structure[wave1data$PFATHER==1 & wave1data$PLIVE==2 & wave1data$PPARTPAR==1]="Two Parent Biological"
# wave1data$family.structure[wave1data$PMOTHER==1 & (wave1data$PSPRELAT==2 | wave1data$PSPRELAT==3)]="Two Parent Nonbiological"
# wave1data$family.structure[wave1data$PFATHER==1 & (wave1data$PSPRELAT==2 | wave1data$PSPRELAT==3)]="Two Parent Nonbiological"
# wave1data$family.structure[(wave1data$PMOTHER==2 | wave1data$PMOTHER==3) & (wave1data$PSPRELAT==1 | wave1data$PSPRELAT==2 | wave1data$PSPRELAT==3)]="Two Parent Nonbiological"
# wave1data$family.structure[(wave1data$PFATHER==2 | wave1data$PFATHER==3) & (wave1data$PSPRELAT==1 | wave1data$PSPRELAT==2 | wave1data$PSPRELAT==3)]="Two Parent Nonbiological"
# wave1data$family.structure[is.na(wave1data$family.structure)]="Single Parent/Other"

wave1data$family.structure[wave1data$PMOTHER==1 & wave1data$PLIVE==1 & wave1data$PSPRELAT==1]=1
wave1data$family.structure[wave1data$PMOTHER==1 & wave1data$PLIVE==2 & wave1data$PPARTPAR==1]=1
wave1data$family.structure[wave1data$PFATHER==1 & wave1data$PLIVE==1 & wave1data$PSPRELAT==1]=1
wave1data$family.structure[wave1data$PFATHER==1 & wave1data$PLIVE==2 & wave1data$PPARTPAR==1]=1
wave1data$family.structure[wave1data$PMOTHER==1 & (wave1data$PSPRELAT==2 | wave1data$PSPRELAT==3)]=1
wave1data$family.structure[wave1data$PFATHER==1 & (wave1data$PSPRELAT==2 | wave1data$PSPRELAT==3)]=1
wave1data$family.structure[(wave1data$PMOTHER==2 | wave1data$PMOTHER==3) & (wave1data$PSPRELAT==1 | wave1data$PSPRELAT==2 | wave1data$PSPRELAT==3)]=1
wave1data$family.structure[(wave1data$PFATHER==2 | wave1data$PFATHER==3) & (wave1data$PSPRELAT==1 | wave1data$PSPRELAT==2 | wave1data$PSPRELAT==3)]=1
wave1data$family.structure[is.na(wave1data$family.structure)]=0

# Highest parent education in household
dadeductemp=rep(NA,nrow(wave1data))
dadeductemp[wave1data$PDADEDUC==0 | wave1data$PDADEDUC==1 | wave1data$PDADEDUC==2]=0
dadeductemp[wave1data$PDADEDUC==3 | wave1data$PDADEDUC==4 | wave1data$PDADEDUC==5 | wave1data$PDADEDUC==7]=1
dadeductemp[wave1data$PDADEDUC==6 | wave1data$PDADEDUC==8]=2
dadeductemp[wave1data$PDADEDUC==9 | wave1data$PDADEDUC==10]=3
dadeductemp[wave1data$PDADEDUC>=11 & wave1data$PDADEDUC<=14]=4
momeductemp=rep(NA,nrow(wave1data))
momeductemp[wave1data$PMOMEDUC==0 | wave1data$PMOMEDUC==1 | wave1data$PMOMEDUC==2]=0
momeductemp[wave1data$PMOMEDUC==3 | wave1data$PMOMEDUC==4 | wave1data$PMOMEDUC==5 | wave1data$PMOMEDUC==7]=1
momeductemp[wave1data$PMOMEDUC==6 | wave1data$PMOMEDUC==8]=2
momeductemp[wave1data$PMOMEDUC==9 | wave1data$PMOMEDUC==10]=3
momeductemp[wave1data$PMOMEDUC>=11 & wave1data$PMOMEDUC<=14]=4
parents.highest.educ=pmax(dadeductemp,momeductemp,na.rm=TRUE)
# wave1data$highest.education.parent.in.household=rep(NA,nrow(wave1data))
# wave1data$highest.education.parent.in.household[parents.highest.educ==0]="Less than high school"
# wave1data$highest.education.parent.in.household[parents.highest.educ==1]="High school degree"
# wave1data$highest.education.parent.in.household[parents.highest.educ==2]="AA/vocational degree"
# wave1data$highest.education.parent.in.household[parents.highest.educ==3]="BA/BS degree"
# wave1data$highest.education.parent.in.household[parents.highest.educ==4]="Higher degree"
# wave1data$highest.education.parent.in.household[is.na(parents.highest.educ)]="Missing"

wave1data$highest.education.parent.in.household=rep(NA,nrow(wave1data))
wave1data$highest.education.parent.in.household[parents.highest.educ==0]=0
wave1data$highest.education.parent.in.household[parents.highest.educ==1]=1
wave1data$highest.education.parent.in.household[parents.highest.educ==2]=1
wave1data$highest.education.parent.in.household[parents.highest.educ==3]=2
wave1data$highest.education.parent.in.household[parents.highest.educ==4]=3
#wave1data$highest.education.parent.in.household[is.na(parents.highest.educ)]=mean(parents.highest.educ,na=T)
wave1data$highest.education.parent.in.household.mis=is.na(parents.highest.educ)

# Gender of teenager
wave1data$gender=rep(NA,nrow(wave1data))
#wave1data$gender[wave1data$TEENSEX==0]="MALE"
#wave1data$gender[wave1data$TEENSEX==1]="FEMALE"
wave1data$female=wave1data$TEENSEX

# Race/ethnicity of teenager
wave1data$race.ethnicity=rep(NA,nrow(wave1data))
# wave1data$race.ethnicity[wave1data$TEENRACE==1]="White/Other"
# wave1data$race.ethnicity[wave1data$TEENRACE==2]="African American"
# wave1data$race.ethnicity[wave1data$TEENRACE==3]="Hispanic"
# wave1data$race.ethnicity[wave1data$TEENRACE>=4]="White/Other"

wave1data$race.black=wave1data$TEENRACE==2
wave1data$race.hispanic=wave1data$TEENRACE==3

# Age of teenager
wave1data$age.teenager=wave1data$AGE
wave1data$age.missing=(wave1data$AGE==888)
# Fill in mean value for teenager with missing age
wave1data$age.teenager[wave1data$AGE==888]=NA
#wave1data$age.teenager[is.na(wave1data$age.teenager)]=mean(wave1data$age.teenager,na.rm=TRUE)

# Has student dropped out of school
wave1data$school.dropout=(wave1data$PSCHTYP==4)

# Work intensity (intensity of adolescent employment)
wave1data$work.intensity=rep(NA,nrow(wave1data))
wave1data$work.intensity[wave1data$WORKHRS==0]="Nonworker"
# Intense: >=20 hours
wave1data$work.intensity[wave1data$WORKHRS>=1 & wave1data$WORKHRS<20]="Moderate"
wave1data$work.intensity[wave1data$WORKHRS>=20 & wave1data$WORKHRS<200]="Intense"

# Alcohol use
wave1data$alcohol.use=rep(NA,nrow(wave1data))
wave1data$alcohol.use[wave1data$DRINK==7]="Never"
wave1data$alcohol.use[wave1data$DRINK==5 | wave1data$DRINK==6]="Moderate"
wave1data$alcohol.use[wave1data$DRINK<=4]="Heavy"

# Marijuana use
wave1data$marijuana.use=rep(NA,nrow(wave1data))
wave1data$marijuana.use[wave1data$POT==1]="Never"
wave1data$marijuana.use[wave1data$POT==2]="Experimenter"
wave1data$marijuana.use[wave1data$POT==3 | wave1data$POT==4]="Continual User"

## Drop from consideration for matching fifth and sixth graders; students missing work intnsity, alcohol use and marijuana use; students with moderate working intensity
wave1data$not.included.in.sample=(wave1data$PSCHGRA2==5 | wave1data$PSCHGRA2==6 | wave1data$age.missing==TRUE | is.na(wave1data$work.intensity) | is.na(wave1data$alcohol.use) | is.na(wave1data$marijuana.use) | wave1data$work.intensity=="Moderate")

# Create variable which identifies whether wave 1 interview exists for subject
interviewerdata=read.csv("C:/Users/ruoqi/Desktop/Penn/research/Dylan-ThickDescription/ivlink.csv")
wave1interviews=interviewerdata$ids[!(interviewerdata$iver=="W3" | interviewerdata$iver=="W4")]
wave1data$wave1.interview=wave1data$IDS%in%wave1interviews

wave1data$wave1.interview=wave1data$wave1.interview& (!wave1data$family.income.mis) & (!wave1data$highest.education.parent.in.household.mis)

data=wave1data
dsub=data[which(data$not.included.in.sample==FALSE),]
dim(dsub) #2816  932
dsub=dsub[which(dsub$work.intensity!='Moderate'),]
dim(dsub) # 2816  932
dsub$intense=rep(0,dim(dsub)[1])
dsub$intense[which(dsub$work.intensity=='Intense')]=1


#propensity score
dsub$family.income.impute=dsub$family.income
dsub$family.income.impute[dsub$family.income.mis==1]=mean(dsub$family.income,na.rm=TRUE)
dsub$highest.education.parent.in.household.impute=dsub$highest.education.parent.in.household
dsub$highest.education.parent.in.household.impute[dsub$highest.education.parent.in.household.mis==1]=mean(dsub$highest.education.parent.in.household,na.rm=T)
model<-glm(intense~family.income.impute+family.income.mis+
             highest.education.parent.in.household.impute+highest.education.parent.in.household.mis+
             female+race.black+race.hispanic+age.teenager+school.dropout,
           family=binomial(link='logit'),data=dsub,x=TRUE)

x=subset(dsub[c('family.income.impute','family.income.mis','family.structure',
                'highest.education.parent.in.household.impute','highest.education.parent.in.household.mis',
                'female','race.black','race.hispanic','age.teenager','school.dropout')])
pred <- predict(model, newdata = x, type = 'response')
dsub$p=pred
dsub$plogit=car::logit(pred)
#boxplot(prop~intense,data=dsub)
dsub=subset(dsub[c('IDS','intense','family.income','family.income.impute','family.income.mis','family.structure',
                'highest.education.parent.in.household','highest.education.parent.in.household.impute','highest.education.parent.in.household.mis',
                'female','race.black','race.hispanic','age.teenager','school.dropout','alcohol.use','marijuana.use','p','plogit')])
nysr=dsub
save(nysr, file = "nysr.rda")
}
\source{
The National Survey of Youth and Religion.
}
\references{
Longest, K. C. and Shanahan M. J., Adolescent Work Intensity and Substance Use: The Mediational and Moderational Roles of Parenting, Journal of Marriage and Family, Vol. 69, No. 3, pp. 703-720.
}
\examples{
data("nysr")
summary(nysr)
}
\keyword{datasets}
