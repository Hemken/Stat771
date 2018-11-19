Admit <- UCBAdmissions

# Overall admissions
apply(Admit, 1, sum)
prop.table(apply(Admit, 1, sum))

# Gender distribution
apply(Admit, 2, sum)
prop.table(apply(Admit, 2, sum))

# Admission ~ Gender
apply(Admit, c(2,1), sum)
prop.table(apply(Admit, c(2,1), sum), margin=1)
chisq.test(apply(Admit, c(2,1), sum))

# Dept distribution
apply(Admit, 3, sum)
prop.table(apply(Admit, 3, sum))

# Admission ~ Gender
apply(Admit, c(3,1), sum)
prop.table(apply(Admit, c(3,1), sum), margin=1)
chisq.test(apply(Admit, c(3,1), sum))
