rm(list = ls())
graphics.off()
#A mixed effects model for a 2AFC paradigm
#"Chicken colour discrimination depends on background colour"
#Peter Olsson, Robin D. Johnsson, James Foster, John Kirwan, Olle Lind and Almut Kelber
#################################################################
#	Useful Functions & Settings									#
#################################################################
#script of my most used functions
source(paste0(Sys.getenv('HOME'),'/Dropbox/R scripts/', 'TMUF.R') )

Instalload(c('rstan','brms'))
#make sure your R version is up to date before trying to load brms
rstan_options(auto_write = TRUE)#; options(mc.cores = parallel::detectCores())
options(mc.cores = parallel::detectCores() - 1)#the model will be estimated on all CPUs except one, so that the user can do other things while they wait

#################################################################
#	Organise the Dataset										#
#################################################################

#Load the data
#Peter's original "short format" summary has total successes in each row
load(file = paste0(Sys.getenv('HOME'),'/Dropbox/Colour discrimination on chromatic backgrounds/Data/for James/', 'colourdata.shortformat-rebatch','.Rdata'))
#My extended "long format" table has a row for each trial with success or failure
load(file = paste0(Sys.getenv('HOME'),'/Dropbox/Colour discrimination on chromatic backgrounds/Data/for James/', 'colourdata.longformat-rebatch','.Rdata'))

#Peter's data, a summary of each individual's performance
head(col.dta)

#I converted it to "long" format. Success in each trial was either TRUE or FALSE
head(cd.long)

#add a row that states whether target was the same or different colour type as the background (e.g. green on green = 'same')
cd.long$target <- ifelse(cd.long$target.same, 'same', 'diff')
cd.long$target <- as.factor(cd.long$target)#this is a factor
cd.long$chick <- as.factor(cd.long$chick)#chick is also a factor

#set reference conditions
levels(cd.long$target)#different is the default reference level
#set same as the reference level
cd.long$target <- relevel(cd.long$target,'same')
levels(cd.long$background)
#keep background_green as reference level

# I also split conditions into two factors with an interaction:
# backgrounds could be green or orange, and targets could be the same colour type as their background, or a different one.
# Finally, I added "chick",  which accounts for the fact that "ind" names shared across experiments actually belong to different chicks ("A1" is like "Svensson" in the world of chicks).

#####################################################################
#	Fit a Non-linear Model with All Fixed Effects					#
#####################################################################

#	SOME PROBLEMS
#1. BASELINE TOO VARIABLE √ constrained
#2. LAPSE SEEMS TOO LARGE √ Kuss et al. Beta(2,20) prior is fine
#3. FIXED EFFECTS ARE FAR TOO SMALL √ coefficients are now free
#4. LOWER ASYMPTOTE ONLY REACHED BELOW COL_DIFF = 0 X unfixable
#5. RANDOM EFFECTS OF WIDTH AND LAPSE NOT IMPLEMENTED √ both now

#re: 5, would have been nice to use:
#	lapse ~ 1 + (1 | gr(chick, dist = 'beta'))
#I have now specified the link function via:
#	inv_logit(lapse) ~ 1 + (1|chick)
# which is equivalent to:
#	lapse ~ 1 + (1 | gr(chick, dist = 'bernoulli'))
# (which doesn't yet exist)

#re: 4, I think Alcalá-Quintana & García-Pérez explicitly state some properties as restrictions for their model.
# It would have been nice to add this parameter:
#	a = threshold - width/2,
#	a~1,
#	prior(a) = Normal(0,10, lb = 0)

# I chose the "threshold,width" parametrisation, to directly model changes in inflection-point threshold and size of non-asymptotic region across conditions and individuals.
# see https://pdfs.semanticscholar.org/f7c1/f70ecb98614dcb64e7bff6be03ad692d5d99.pdf
# , http://journalofvision.org/5/5/8/
# and https://www.sciencedirect.com/science/article/pii/S0042698918300567 for nice explanations, 
# and see https://link.springer.com/content/pdf/10.3758/BF03194547.pdf for the theoretical problem this addresses.

#####################################################################
#	Make a Model Formula											#
#####################################################################

#formulae and factors are arranged so that reference condition is:
# background_green & target_same (green on green)
#this is determined by both the reference levels in the data frame.
# The order of "background" and "target" in the formula determines naming
# and which one's independent effects are estimated first.
# i.e. background*target: 1st background, 2nd target, 3rd background:target.
# With broad unbiased priors this effect is negligible.

#threshold and width are estimated on a log scale, which keeps them positive and allows free reign for random effects
#lapse is estimated on a logit scale, also allowing free estimation of random effects in [0,1] space

#model 1, intercept is background_green, target_same 
modnm1 <- 'TW.Model-bg_target_REVISION'#distinguish it from others
frm1 <-              bf(
	formula = success ~ base	 +	#guess rate
						(1-inv_logit(lapse)-base) 	* #curve region
			inv_logit(0++4.39*(	Colour.difference-exp(threshold)	) 	/
						(	exp(width)	)), #threshold-width curve
	  base ~ 1, #baseline has a single mean
	  lapse ~ 1 + (1|chick) +(1|batch),  #lapse rate depends on chick
	  #threshold coef depend on fixef & chick
	 threshold ~ background*target*sex +(1|chick) +(1|batch),
	 #width coef depend on fixef & chick
	 width ~ background*target*sex +(1|chick) +(1|batch),
      nl = TRUE)#the joint distribution for these parameters is undefined, and therefore the parameters themselves are "nonlinear"
      
#####################################################################
#	Select Some Priors												#
#####################################################################

#	Model1, Raneff of lapse; Intercept:target_green&background_green	#

#remember, we need
get_prior(frm1, data = cd.long)[,c('class', 'coef','nlpar', 'group')]#print priors to comm window as a reminder
# these ones are well behaved
# only base has a narrow prior
prr1 <- c(
  #very restrictive prior for guess rate, centred on 0.5
	prior(beta(250,250), nlpar= 'base', lb = 0.25, ub = 0.75),
  #lapse rate is unbiased, but cannot be more than 27%
	prior(normal(-3,10), nlpar= 'lapse', ub = -1),
  # use the default prior for random effects of lapse:
  #	student_t(3, 0, 10)),
  # this can be done by leaving:
  #	prior(... nlpar= 'lapse', class= sd), unassigned
  #Both threshold and width should be positive numbers, probably ≈1
  #i.e. exp(0) = 1
  #beware, bounds on threshold and width priors
  #affect their coefficients (so don't apply bounds)
	prior(normal(0,10), nlpar= 'threshold', class = 'b'),
	prior(normal(0,10), nlpar = 'width', class = 'b'),
  #Coefficient parameters, centred on 0
  #(<0 = param smaller, >0 = param larger)
	prior(normal(0,10), nlpar= 'threshold', coef= 'backgroundorange'),
	prior(normal(0,10), nlpar= 'threshold', coef= 'targetdiff'),
	prior(normal(0,10), nlpar= 'threshold', coef= 'backgroundorange:targetdiff'),
  # use the default prior for random effects of threshold, unassigned
	prior(normal(0,10), nlpar= 'width', coef= 'backgroundorange'),
	prior(normal(0,10), nlpar= 'width', coef= 'targetdiff'),
	prior(normal(0,10), nlpar= 'width', coef= 'backgroundorange:targetdiff')#,
  # use the default prior for random effects of width, unassigned
)

#####################################################################
#	Inspect Stan Code												#
#####################################################################
stc1 <- make_stancode( 
			 formula = frm1,
 	         data = cd.long, family = bernoulli("identity"), 
             prior = prr1
             )#seems to work

stc1#I don't quite know what all of this means
write.table(stc1, file = paste0(Sys.getenv('HOME'),'/Dropbox/Colour discrimination on chromatic backgrounds/Manuscript/Supplemental files/', modnm1,'.stan'), col.names = F, row.names = F, quote = F)			

#####################################################################
#	Run the Model													#
#####################################################################

TW_col.dta1 <- brm( formula = frm1,
 	         data = cd.long, family = bernoulli("identity"), 
             prior = prr1,
             #finely sampled
        		control = list(adapt_delta = 0.99),
        		#IMPORTANT,
        		#random effects of threshold and width
        		#should change together (ideally, explicit correllation).
        		#setting initial values to 0 is a work-around
        		inits = 0,
        		#200 or 2000 iterations give similar results to 10000
        		iter = 1000#10000
                          )
                          
#make sure to save the model object, it is quite large
 save(TW_col.dta1, file = paste0(Sys.getenv('HOME'),'/Dropbox/Colour discrimination on chromatic backgrounds/Manuscript/Supplemental files/', modnm1,'.Rdata'))

#sometimes I automatically quit for the sake of my poor computer
# q('no')

#here's one I made earlier
load(file = paste0(Sys.getenv('HOME'),'/Dropbox/Colour discrimination on chromatic backgrounds/Manuscript/Supplemental files/', modnm1,'.Rdata'))

#inspect the model
summary(TW_col.dta1)#;dev.new(height = 5, width = 7);plot(TW_col.dta1)
marginal_effects(TW_col.dta1)

#to plot, use "Model Plotter 20190603 bg_target REVISION.R"