hist(exposure_outcome$solar$pfoa)

nrow(exposure_outcome$solar)

exposure_outcome$solar <- exposure_outcome$solar |>
  mutate(ses_num = ses |> 
           as.factor() |>
           as.numeric(), 
         ses_num = if_else(ses_num == 4, 0, ses_num),
         ses_num = if_else(ses_num == 5, 4, ses_num))
 reduced <- exposure_outcome$solar |> 
   filter(ses_num != 4, !is.na(pfos))

table(exposure_outcome$solar$ses, 
      exposure_outcome$solar$ses_num)

car::Anova(lm(log(pfos)  ~ as.factor(ses_num) + sex, data = exposure_outcome$solar))
car::Anova(lm(log(pfhxs) ~ as.factor(ses_num) + sex, data = exposure_outcome$solar))
car::Anova(lm(log(pfoa)  ~ as.factor(ses_num) + sex, data = exposure_outcome$solar))
car::Anova(lm(log(pfna)  ~ as.factor(ses_num) + sex, data = exposure_outcome$solar))
car::Anova(lm(log(pfda)  ~ as.factor(ses_num) + sex, data = exposure_outcome$solar))


summary(lm(pfos  ~ as.numeric(ses_num), data = exposure_outcome$solar |> filter(ses_num != 4)))
0.0003137
summary(lm(pfhxs ~ as.numeric(ses_num), data = exposure_outcome$solar |> filter(ses_num != 4)))
0.007086
summary(lm(pfhps ~ as.numeric(ses_num), data = exposure_outcome$solar |> filter(ses_num != 4)))
0.003558
summary(lm(pfoa  ~ as.numeric(ses_num), data = exposure_outcome$solar |> filter(ses_num != 4)))
0.00129
summary(lm(pfna  ~ as.numeric(ses_num), data = exposure_outcome$solar |> filter(ses_num != 4)))
0.01731
summary(lm(pfda  ~ as.numeric(ses_num), data = exposure_outcome$solar |> filter(ses_num != 4)))
0.006953


0.0003137
0.006953
0.007086
0.00129
0.01731


cor.test(x = reduced$pfos,  y = reduced$ses_num)
cor.test(x = reduced$pfhxs, y = reduced$ses_num)
cor.test(x = reduced$pfoa,  y = reduced$ses_num)
cor.test(x = reduced$pfna,  y = reduced$ses_num)
cor.test(x = reduced$pfda,  y = reduced$ses_num)


summary(lm(pfhxs ~ as.numeric(ses_num), data = exposure_outcome$solar |> filter(ses_num != 4)))
summary(lm(pfoa  ~ as.numeric(ses_num), data = exposure_outcome$solar |> filter(ses_num != 4)))
summary(lm(pfna  ~ as.numeric(ses_num), data = exposure_outcome$solar |> filter(ses_num != 4)))
summary(lm(pfda  ~ as.numeric(ses_num), data = exposure_outcome$solar |> filter(ses_num != 4)))


data <- exposure_outcome$solar |>
  mutate(ses_missing = ses == "missing")


chisq.test(data$ses_missing, data$sex)
summary(lm(age ~ ses_missing, data))
summary(lm(age ~ ses_missing, data))
summary(lm(bmi ~ ses_missing, data))
t.test(data$age ~ data$ses_missing)
