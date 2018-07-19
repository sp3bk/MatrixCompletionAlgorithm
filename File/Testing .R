ab <- sample(1:10,2)
ba <- sample(1:10,10)

ba <- append(ba, ab)
ba
saveRDS(ba,file="ab.Rda")



  