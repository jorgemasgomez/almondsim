# Integer rounding shared by maternal and paternal OCS roles.
optisel_round_roles <- function(role,crosses=20L) {
 round_role<-function(indices) {
   desired<-crosses*role$contributions[indices]/sum(role$contributions[indices])
   n<-floor(desired);missing<-crosses-sum(n)
   if(missing>0) {chosen<-head(order(-(desired-n),seq_along(n)),missing);n[chosen]<-n[chosen]+1L}
   n
 }
 mother_counts<-round_role(role$mothers);father_counts<-round_role(role$fathers)
 mothers<-role$index[role$mothers];fathers<-role$index[role$fathers]
 if(identical(mothers,fathers)) {
   desired<-2*crosses*(role$contributions[role$mothers]+role$contributions[role$fathers])
   counts<-floor(desired);missing<-2*crosses-sum(counts)
   if(missing>0) {chosen<-head(order(-(desired-counts),seq_along(counts)),missing);counts[chosen]<-counts[chosen]+1L}
   mother_counts<-floor(counts/2)
   odd<-which(counts%%2==1)
   needed<-crosses-sum(mother_counts)
   if(needed>0) {
     order_odd<-odd[order(-(2*crosses*role$contributions[role$mothers][odd]-mother_counts[odd]),odd)]
     mother_counts[head(order_odd,needed)]<-mother_counts[head(order_odd,needed)]+1L
   }
   father_counts<-counts-mother_counts
 }
 list(mothers=mothers,fathers=fathers,mother_counts=mother_counts,father_counts=father_counts)
}
