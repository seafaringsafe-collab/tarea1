pob = 1:6

res = c()
start.time <- Sys.time()

for (j in 1:100){
  
  dados = sample(pob, 2, replace = TRUE)   # lanza 2 dados
  suma = sum(dados)                        # suma de los dos dados
  
  res <- c(res, suma)
  
  plot(0,0,xlim=c(-50,50), ylim=c(-50,50), axes = FALSE, type='n', xlab = "", ylab = "")
  text(0,10, dados[1], cex=10)             # dado 1
  text(0,-10, dados[2], cex=10)            # dado 2
  text(0,30, paste("Suma =", suma), cex=2) # suma
  
  Sys.sleep(1)
}

end.time <- Sys.time()
time.taken <- end.time - start.time
print(time.taken)

hist(res, breaks=2:13, right = F)