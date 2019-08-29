load("/home/lucy/bin/crema//cpat_models/ath_logit.RData")
test <- read.table(file="HL/LL/dn/names.cpat.dat",sep="\t",col.names=c("ID","mRNA","ORF","Fickett","Hexamer"))
test$prob <- predict(mylogit,newdata=test,type="response")
attach(test)
output <- cbind("mRNA_size"=mRNA,"ORF_size"=ORF,"Fickett_score"=Fickett,"Hexamer_score"=Hexamer,"coding_prob"=test$prob)
write.table(output,file="HL/LL/dn/names.cpat",quote=F,sep="\t",row.names=ID)
