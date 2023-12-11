library(tidyverse)
library(meta)

####导入文件（排序有格式）####
calculation<-read_csv('CALCULATION.csv')

####计算随机效应模型####
resultVER_IVW_random <- data.frame(Outcome = character(), OR = numeric(), LB = numeric(), UB = numeric())
j=1
#i的上限不同文件要改
for (i in 2:30) {

  is_different <- !identical(calculation$Outcomes[i], calculation$Outcomes[i - 1])
  if (is_different){

    extract <- calculation[j:(i-1), ]
    #print(extract)
    m.gen_bin <- metagen(TE = log(OR),
                         lower = log(LB),
                         upper = log(UB),
                         data = extract,
                         sm = "OR",
                         method.tau = "PM",
                         fixed = F,
                         random = T,
                         title = as.character(extract[1,1]))
    
    summary(m.gen_bin)
    resultVER_IVW_random <- rbind(resultVER_IVW_random, data.frame(Outome=m.gen_bin[["title"]],
                                                     OR = exp(m.gen_bin[["TE.random"]]), 
                                   LB = exp(m.gen_bin[["lower.random"]]), 
                                   UB = exp(m.gen_bin[["upper.random"]])))
    j=i
  }
}

resultVER_IVW_random$type<-"random"
print(resultVER_IVW_random)
#write.csv(resultVER_IVW_random, file = "reverse_result_WM_random.csv", row.names = T)


####计算固定效应模型####
resultVER_IVW_fix <- data.frame(Outcome = character(), OR = numeric(), LB = numeric(), UB = numeric())
j=1
#i的上限不同文件要改
for (i in 2:30) {
  
  is_different <- !identical(calculation$Outcomes[i], calculation$Outcomes[i - 1])
  if (is_different){
    
    extract <- calculation[j:(i-1), ]
    #print(extract)
    m.gen_bin <- metagen(TE = log(OR),
                         lower = log(LB),
                         upper = log(UB),
                         data = extract,
                         sm = "OR",
                         method.tau = "PM",
                         fixed = T,
                         random = F,
                         title = as.character(extract[1,1]))
    
    summary(m.gen_bin)
    resultVER_IVW_fix <- rbind(resultVER_IVW_fix, data.frame(Outome=m.gen_bin[["title"]],
                                                                   OR = exp(m.gen_bin[["TE.common"]]), 
                                                                   LB = exp(m.gen_bin[["lower.common"]]), 
                                                                   UB = exp(m.gen_bin[["upper.common"]])))
    j=i
  }
}
 
resultVER_IVW_fix$type<-"fixed"
print(resultVER_IVW_fix)
#write.csv(resultVER_IVW_fix, file = "reverse_result_WM_fix.csv", row.names = T)




