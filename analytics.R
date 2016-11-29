library(gdata);
library('gridExtra');
library('grid');
library('gtable');
library('agricolae');
library('multcomp');
library(plotrix);
library('nlme');


# R version 2.14.1 (2011-12-22)
source('M:/anyela/repo/senescence_disease/generic.R');

setwd("M:/anyela/repo/root")

read_file = function(filename)
{
  d = read.table(filename, head=TRUE, sep=',');
  return(d);
  
}


readData = function(phenofile, metafile)
{
  f = list();
  f$p = read_file(phenofile);
  
  f$m = read_file(metafile);

  return(f);
}

# Convert to two cols
createMatrix = function(data, kw, mx)
{
  
  s = data[which(data$S.W == kw),];
  u = c();
  for(i in 1:nrow(s))
  {
    u = rbind(rbind(u, cbind(s[i,1:2], dist=seq(5,(ncol(s)-3)*5,by=5))));
  }
  
  l = unmatrix(as.matrix(s[,4:ncol(s)]), byrow = TRUE);
  if(!is.na(mx)) { l = rescale(unmatrix(as.matrix(s[,4:ncol(s)]), byrow = TRUE), c(1,mx)); }
  
  m = cbind(u,l);
  colnames(m)[ncol(m)] = tolower(kw);
  return(m);
}

createFrame = function(data, mx)
{
  s = createMatrix(data, 'S', mx)
  w = createMatrix(data, 'W', mx);
  m = cbind(s,w);
  m = m[,c('Plant', 'day', 'dist', 's', 'w')];
  return(m);
  
}

create_corrmat = function(data)
{
  y = aggregate(cbind(s, w) ~ Plant + dist, FUN=mean, data=data);
  o = y[order(y$Plant, y$dist), ];
  
  h = list()
  for(pname in unique(o$Plant))
  {
    a = as.character(pname)
    subdata = data[which(o$Plant == pname),c(3:4)];
    l = as.numeric(unmatrix(as.matrix(subdata), byrow = FALSE));
    h[[a]] = l;
  }
  
  return(data.frame(h));
}


calculateCorrelation = function(data, trait)
{
  
  m = matrix(NA, nrow=50, 50);
  for(i in 1:50)
  {
    for(j in 1:50)
    {
      p1 = which(copydata$PopNo == i);
      p2 = which(copydata$PopNo == j);
      m[i,j] = cor(copydata[p1,'s'], copydata[p2,trait])
    }
  }
  colnames(m) = paste('p', 1:50, sep='');
  rownames(m) = paste('p', 1:50, sep='');
  heatmap(m, main = trait);
  #return(m);
}


interactionplots = function(mx)
{
  png('interactionWDist.png', width = 1200, height = 800,res=200);
  with(m2, interaction.plot(PopNo, dist, w, col=rainbow(5), lwd=2));
  dev.off();
  
  png('interactionSDist.png', width = 1200, height = 800,res=200);
  with(m2, interaction.plot(PopNo, dist, s, col=rainbow(5), lwd=2));
  dev.off();
  
  png('interactionWPop.png', width = 1500, height = 1600,res=200)
  with(m2, interaction.plot(dist, PopNo, w, col=rainbow(50), lwd=2, lty=c(1:51),legend=F));
  legend("topright", legend=unique(m2$PopNo),bty="n",lty=c(1:51),lwd=2, 
         col=rainbow(50), cex=0.55,title="Population",inset = 0) 
  dev.off();
  
  png('interactionSPop.png', width = 1500, height = 1600,res=200)
  with(m2, interaction.plot(dist, PopNo, s, col=rainbow(50), lwd=2, lty=c(1:51), legend=F));
  legend("topright", legend=unique(m2$PopNo),bty="n",lty=c(1:51),lwd=2, 
         col=rainbow(50), cex=0.55,title="Population",inset = 0) 
  dev.off();
  copydata = aggregate(cbind(s,w) ~  dist + PopNo, data=m2, 
                       FUN= "mean", na.action = na.pass);
  
  png('HMS.png', width = 1000, height = 1000,res=150);
  calculateCorrelation(copydata, 's')
  dev.off();
  
  png('HMW.png', width = 1000, height = 1000,res=150);
  calculateCorrelation(copydata, 'w')
  dev.off();
  
  m = aov(s ~ genotype, m2);
  o = with(m2, pairwise.t.test(s, genotype, p.adj = "bonf"));
  stable = data.frame(as.matrix(o$p.value))
  
  m = aov(w ~ genotype, m2);
  o = with(m2, pairwise.t.test(w, genotype, p.adj = "bonf"));
  wtable = data.frame(as.matrix(o$p.value))
  
 
  
}

plotInteraction = function(e, ct, ge, gc, data, dt, alttitle=NULL)
{
  print(dt)
  copydata = data;
  cr = length(unique(copydata[[dt]]));
  for(experimentline in e)
  {
    i = which(copydata[[ge]] == experimentline)
    sexp  = copydata[i,];
    for(controlline in ct)
    {
      par(mfrow=c(2,2))
      i = which(copydata[[gc]] == controlline)
      scont  = copydata[i,];
      mx = copydata[union(which(copydata[[ge]] == experimentline), which(copydata[[ge]] == controlline)),];
      for(trait in c('s', 'w'))
      {
        x = max(mx[[trait]]); y = min(mx[[trait]])
        interaction.plot(sexp[['date']], sexp[[dt]], sexp[[trait]], col=rainbow(cr), 
                         main = paste(experimentline, alttitle, sep=' '), ylab= trait, xlab='month', 
                         lwd=2, legend=F, cex.axis=1, ylim=c(y, x-1));
        legend("topleft", legend=unique(sexp[[dt]]), bty="n",lwd=2, 
               col=rainbow(cr), cex=0.80,title="dist",inset = 0);
        
        interaction.plot(scont[['date']], scont[[dt]], scont[[trait]], col=rainbow(cr), 
                         main = controlline, ylab= trait, xlab='month', lwd=2, ylim=c(y, x-1), legend=F);
        legend("topleft", legend=unique(scont[[dt]]), bty="n",lwd=2, 
               col=rainbow(cr), cex=0.80,title=dt,inset = 0);
        
      }
    }
  }
}

 
postHoc = function(e=t1exp1, ct=t1control1, ge='genotype', gc='genotype', data=m, trait, trait2, tp,
                   om, ma)
{
  colbox=c("pink", "red", "blue", 'green', 'black', 'brown', 'yellow', 'orange');
  l = c();
  for(experimentline in e)
  {
    for(controlline in ct)
    {
      par(mfrow=c(3,2), oma=om, mar=ma)
      for(dt in unique(data[[trait]]))
      {
        i = intersect(union(which(data[[ge]] == experimentline), 
                               which(data[[gc]] == controlline)),
                               which(data[[trait]] == dt));
        sub = data[i,];
        intmodel = interaction(sub[[trait2]], sub[[ge]]);
        y = as.formula(paste(tp, '~ intmodel', sep =''))
        m = aov(y, sub);
        tuk <- glht(m, linfct = mcp(intmodel = "Tukey"));
        t = summary(tuk);
        tuk.cld <- cld(tuk)   # letter-based display
        j = length(unique(unique(data[[trait2]])));
        plot(tuk.cld, col=colbox[1:j], cex.axis=0.5, las=2, xaxt="n", 
             xlab="", ylab= paste(tp, 'at ', dt), ylim=c(0, max(sub[[tp]])+1));
        
        text(1.5,0.5, cex=1, controlline); text(j+1.5,0.5, experimentline, cex=1);
        l <- rbind(l, cbind(dt, as.matrix(t$test$tstat), p.adj=t$test$pvalues));
        abline(v=j+0.5)
        axis(1, at=c(1:(j*2)),labels=rep(unique(data[[trait2]]),2), cex.axis=0.7);
      }
      #title(paste(trait, experimentline, 'vs', controlline), outer=TRUE)
    }
  }
  
  csep=' - ';
  msep='\\.'
  l = data.frame(test=rownames(l),l, row.names = 1:nrow(l) );
  l = data.frame(l, monthA=sapply(l$test, function(x) processRow(x,1,1,csep, msep)));
  l = data.frame(l, GenotypeA=sapply(l$test, function(x) processRow(x,1,2,csep, msep)));
  l = data.frame(l, monthB=sapply(l$test, function(x) processRow(x,2,1,csep, msep)));
  l = data.frame(l, GenotypeB=sapply(l$test, function(x) processRow(x,2,2,csep, msep)));
  l = data.frame(l, samemonth=sapply(rownames(l), 
                                     function(x) compareValues(as.character(l$monthA[as.numeric(x)]), 
                                                               as.character(l$monthB[as.numeric(x)])) ));
  return(l)
}



# split rowname
# n1 = 1 indicates left of '-' and 2 right
# n2 = 1 indicates lrft of '.' and 2 right
processRow = function(x, n1, n3, csep, msep)
{
  x= as.character(x);
  s = strsplit(strsplit(x, csep)[[1]][n1], msep)[[1]][n3]
  return(s[[1]])
  
}


testdummy = function()
{
  i = which(m2$genotype == 'subline1');
  sexp  = m2[i,];
  j = which(m2$genotype == 'Meltra');
  scont  = m2[j,];
  trait='s'
  ge = 'genotype'
  gc = 'genotype'
}

changeFieldName = function(data, lfname)
{
  
  data = data.frame(data, idtag = sapply(data$Plant, function(x) 
    if(x < 10000) { paste(lfname, '-','0',x, sep='')} 
    else { paste(lfname, '-',x, sep='')}))
  return(data);
}

# Create name based on metadata
createID = function(data, metadata, rawid, lfname)
{
  s = as.numeric(strsplit(as.character(rawid/10), '\\.')[[1]][1]);
  i = which(metadata[['NPPC.Plant.no']] == s[[1]][1])[1];

  if(s < 10) 
  { 
    n = paste(lfname, '-','0',s, 
          metadata[['roots']][i],  metadata[['treatment']][i], metadata[['Rep.1']][i], sep='')
  } 
  else 
  { 
    n = paste(lfname, '-', s, 
          metadata[['roots']][i], metadata[['treatment']][i], metadata[['Rep.1']][i], sep='')
  }
  return(n)
  
}


changeFieldNameLF = function(data, metadata, lfname)
{
  
  data = data.frame(data, idtag = sapply(data$Plant, function(x) createID(data, metadata,x, lfname)));
  return(data);
}

# change a string fro another
changeWord = function(kw, word_list)
{
  #print(kw)
  return(which(names(word_list) == kw))
}

changeDate = function(data, date_list)
{
  data = data.frame(data, date=sapply(data[['day']], function(x) changeWord(x,date_list)))
  return(data);
}


processPhaseOne = function()
{
  g = readData('LF5_pheno.csv', 'LF5P1_meta.csv');
  m = createFrame(g$p);
  f = cor(create_corrmat(m));
  #heatmap(f)
  
  m = changeFieldName(m, 'LF5');
  m = merge(m, g$m, by.x='idtag', by.y='barcode');
  m = changeDate(m, list('351'=1, '379' =2, '408'=3, '450'=4));
  
  colnames(m)[21] = 'genotype'
  m = m[order(m$idtag, m$Plant, m$date, m$dist),];
  m = convert2factor(m, c('date', 'dist', 'Plant', 'rep', 'ColPos', 'PopNo'));
  
  t1exp1 = as.character(unique(m$genotype)[9:15]);
  t1control1 = as.character(unique(m$genotype)[5:8]);
  t1exp2 = as.character(unique(m$PopNo)[34:51]);
  t1control2 = as.character(unique(m$genotype)[c(5:8,11)]);
  t1exp3 = as.character(unique(m$genotype)[1:4]);
  t1control3 = as.character(unique(m$genotype)[6:7]);
  
  break()
  pdf('T1intplot.pdf')
  plotInteraction(t1exp1, t1control1, 'genotype', 'genotype', dist,  m);
  plotInteraction(t1exp2, t1control2, 'PopNo', 'genotype', dist, m, 'recomb_subline3');
  plotInteraction(t1exp3, t1control3, 'genotype', 'genotype', dist, m);
  dev.off()
  
  pdf('T1posthoc1_s.pdf')
  p1=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'w');
  dev.off();
  pdf('T1posthoc2_s.pdf')
  p2=postHoc(t1exp2, t1control2, 'PopNo', 'genotype', m, 'recomb_subline3');
  dev.off();
  pdf('T1posthoc3_s.pdf')
  p3=postHoc(t1exp3, t1control3, 'genotype', 'genotype', m);
  dev.off();
  
  pf = rbind(p1, p2, p3)       
  
  write.table(stable[,-1], file = 'resultsint_s.csv', sep=',');
  write.table(wtable[,-1], file = 'resultsint_w.csv', sep=',');
  write.table(pf, file='pairwise_s.csv', sep=',', quote = F, row.names = F);
  
  fm = aov(s ~ day + dist + genotype, m);
  summary(fm);
  plot(fm);
  TukeyHSD(fm, 'dist');
  TukeyHSD(fm, 'day');
  TukeyHSD(fm, 'genotype')
  
  fm1 = aov(s ~ tx, m);
  r = HSD.test(fm1, "tx", group=TRUE)
  
  summary(tuk)          # standard display
}

checkOutlier = function(data, item_list)
{
  ol = c();
  for(t in c('s', 'w'))
  {
    for(genotypename in item_list)
    {
      i = which(data[['genotype']] == genotypename)
      sub = data[i,];
      #boxplot(sub[['s']] ~ sub[['date']], main=genotypename );
      outliers = boxplot(sub[[t]] ~ sub[['date']], plot=FALSE)$out
      ol = append(ol,rownames(sub[sub[[t]] %in% outliers,]) );
    }
  }
  ol = as.numeric(ol)
  data = data[-ol, ];
  rownames(data) = 1:nrow(data);
  return(data)
}


compareValues = function(A,B)
{
  
  if(A == B)
    return('1')
  else
    return('0');
}

findH = function(data, klist)
{
  g1 = intersect(which(data[['genotype']] == as.character(klist[1]$GenotypeA)), 
                 which(data[['date']] == as.numeric(klist[2]$monthA)));
  g2 = intersect(which(data[['genotype']] == as.character(klist[3]$GenotypeB)), 
                 which(data[['date']] == as.numeric(klist[4]$monthB)));
  mean(data[g1, 's'])
  mean(data[g2, 's'])
  if(mean(data[g1, 's']) >  mean(data[g2, 's']))
  {
    return(1)
  }
  if(mean(data[g1, 's']) <  mean(data[g2, 's']))
  {
    return(2)
  }
  else
    return(0)
}

extractPvalue = function(l, data, csep, msep)
{
  
  
  p = data.frame(l$test$pvalues);
  p = data.frame(test=rownames(p),p, row.names = 1:nrow(p) );
  p = data.frame(p, GenotypeA =sapply(p$test, function(x) processRow(x,1,1,csep, msep)));
  p = data.frame(p, monthA=sapply(p$test, function(x) processRow(x,1,2, csep, msep)));
  p = data.frame(p, GenotypeB =sapply(p$test, function(x) processRow(x,2,1,csep, msep)));
  p = data.frame(p, monthB=sapply(p$test, function(x) processRow(x,2,2, csep, msep)));
  p = data.frame(p, dist1=sapply(p$test, function(x) processRow(x,1,3, csep, msep)));
  p = data.frame(p, dist2=sapply(p$test, function(x) processRow(x,2,3, csep, msep)));
  
  p = data.frame(p, row.names=1:nrow(p));
  p = data.frame(p, samemonth=sapply(rownames(p), 
                                     function(x) compareValues(as.character(p$monthA[as.numeric(x)]), 
                                                               as.character(p$monthB[as.numeric(x)])) ));
  p = data.frame(p, samedist=sapply(rownames(p), 
                                     function(x) compareValues(as.character(p$dist1[as.numeric(x)]), 
                                                               as.character(p$dist2[as.numeric(x)])) ));
  
  p1 = data.frame(p, h=sapply(as.numeric(rownames(p)), function(x) 
    findH(data, p[x, c("GenotypeA","monthA","GenotypeB","monthB")])));
  h=sapply(as.numeric(rownames(p)[1]), function(x) 
    findH(data, p[x, c("GenotypeA","monthA","GenotypeB","monthB")]));
  
  return(l);
}


fitModel = function(data)
{
  data = m;
  data$tall = with(data, interaction(genotype, date,subdist, drop=T, sep = "_"))
  fm = lme(w ~ tall, data=data, random = ~1|Plant, method='ML');
  anova(fm)
  l = summary(glht(fm, linfct=mcp(tall = "Tukey")), test = adjusted(type = "bonferroni"))
  csep=' - ';
  msep='_'
  ll = extractPvalue(data, l, csep, msep);
  write.table(p, file='testsubdist.csv', sep=',');
}

hsdPairWise = function()
{
  library(foreign)
  library(multcomp)
  yield <- read.dta("http://www.stata-press.com/data/r12/yield.dta")
  tx <- with(yield, interaction(fertilizer, irrigation))
  amod <- aov(yield ~ tx, data=yield)
  tuk <- glht(amod, linfct = mcp(tx = "Tukey"))
  
  summary(tuk)          # standard display
  tuk.cld <- cld(tuk)   # letter-based display
  opar <- par(mai=c(1,1,1.5,1))
  plot(tuk.cld)
  par(opar)
  library(agricolae)
  HSD.test(amod, "tx", group=TRUE)
  
}

createSubGroup = function(data)
{
  data = data.frame(data, subdist = sapply(data[['dist']], 
                                   function(x) 
                                     if(as.integer(as.character(x)) >= 20) { 1 } else {0}));
  return(data);
}


processPhaseTwo = function()
{
  g = readData('LF7_pheno.csv', 'LF7P2a_meta.csv');
  m = createFrame(g$p, mx=NA);
  m = changeFieldNameLF(m, g$m, 'LF7');
  m = merge(m, g$m, by.x='idtag', by.y='barcode');
  m = changeDate(m, list('14/04/2015'= 1, '12/05/2015' = 2, '09/06/2015'= 3, '16/07/2015'= 4));
  
  colnames(m)[21] = 'genotype'
  m = m[order(m$idtag, m$Plant, m$date, m$dist),];
  m = convert2factor(m, c('dist', 'Rep', 'ColPos', 'PopNo'));
  #m = checkOutlier(m, unique(m$genotype));
  m = m[c(-506, -538), ];
  
  o = which(m$genotype == 0)
  m = m[-o,];
  m = resetlevels(m, 'genotype')
  m = createSubGroup(m);
  return(m)
}
 



plotPhaseTwo = function()
{
  t1exp1 = c('Bx510', 'Bx509');
  t1control1 = c('AberStar');
  pdf('T1intplot.pdf')
  plotInteraction(t1exp1, t1control1, 'genotype', 'genotype',  m, 'dist',);
  plotInteraction(t1exp1, t1control1, 'genotype', 'genotype', m, 'subdist');
  dev.off()
  
  pdf('T1posthoc1_s_dist.pdf')
  pdist1=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'dist', 'date', 's', c(0,0,2,0), c(4,4,6,1));
  dev.off();
  
  pdf('T1posthoc1_s_subdist.pdf')
  pdist1=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'subdist', 'date', 
                 's', c(0,0,2,0), c(4,4,6,1));
  dev.off();
  write.table(pdist1, file='resPHsubdist_g1.csv', sep=',')

  ##
  pdf('T1posthoc1_s_date.pdf')
  pdate1=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'date', 'dist', 's', c(0,0,2,0), c(4,4,6,1));
  dev.off();
  
  ##
  
  pdf('T1posthoc1_w_dist.pdf')
  pdist1w=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'dist', 'date', 'w', c(0,0,2,0), c(4,4,6,1));
  dev.off();
  
  pdf('T1posthoc1_w_date.pdf')
  pdate1w=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'date', 'dist', 'w', c(0,0,2,0), c(4,4,6,1));
  dev.off();
  
  t1exp1 = c('Bx511', 'Bx514');
  t1control1 = c('AberBite');
  pdf('T2intplot.pdf')
  plotInteraction(t1exp1, t1control1, 'genotype', 'genotype', m, 'dist');
  dev.off()
  
  pdf('T2posthoc1_s_dist.pdf')
  pdist2=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'dist', 'date', 's', c(0,0,2,0), c(4,4,6,1));
  dev.off();
  
  pdf('T2posthoc1_s_subdist.pdf')
  pdist2=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'subdist', 'date', 's', c(0,0,2,0), c(4,4,6,1));
  dev.off();
  write.table(pdist2, file='resPHsubdist_g2.csv', sep=',')
    
  pdf('T2posthoc1_s_date.pdf')
  pdate2=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'date', 'dist', 's', c(0,0,2,0), c(4,4,6,1));
  dev.off();
 
  
  pdf('T2posthoc1_w_dist.pdf')
  pdist2w=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'dist', 'date', 'w', c(0,0,2,0), c(4,4,6,1));
  dev.off();
  
  pdf('T2posthoc1_w_date.pdf')
  pdate2w=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'date', 'dist', 'w', c(0,0,2,0), c(4,4,6,1));
  dev.off();
  
  o = which(m$genotype == 0)
  m = m[-o,];
  m = resetlevels(m, 'genotype')
  pdf('cor_s_w_date.pdf')
  ggplot(m, aes(x=w, y=s, colour=date))  + geom_point() + facet_grid(. ~ genotype)
  dev.off()
  pdf('cor_s_w_dist.pdf')
  ggplot(m, aes(x=w, y=s, colour=dist))  + geom_point() + facet_grid(. ~ genotype)
  dev.off()
  
  t1exp1 = c('Bx510', 'Bx514');
  t1control1 = c('Bx509');
  pdf('T3posthoc1_s_date.pdf')
  pdate1=postHoc(t1exp1, t1control1, 'genotype', 'genotype', m, 'date', 'dist', 's', c(0,0,2,0), c(4,4,6,1));
  dev.off()
  
}


m = processPhaseTwo();
