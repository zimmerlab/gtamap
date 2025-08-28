package server

type GeneralInfo struct {
	NumReadsTotal      int `json:"numReadsTotal"`
	NumReadsAccepted   int `json:"numReadsAccepted"`
	NumReadsDiscarded  int `json:"numReadsDiscarded"`
	NumReadsInProgress int `json:"numReadsInProgress"`
	NumReadsTodo       int `json:"numReadsTodo"`
}

func (h *MappingDataHandler) GetGeneralInfo() *GeneralInfo {

	accepted := 0
	inProgress := 0
	discarded := 0
	for _, qclust := range h.QnameCluster {
		if qclust.ClusterR1.AcceptedRecord != nil || qclust.ClusterR2.AcceptedRecord != nil {
			if qclust.ClusterR1.AcceptedRecord != nil && qclust.ClusterR2.AcceptedRecord != nil {
				accepted++
			} else {
				inProgress++
			}
		}
		if qclust.IsDiscarded {
			discarded++
		}
	}

	g := &GeneralInfo{
		NumReadsTotal:      len(h.QnameCluster),
		NumReadsAccepted:   accepted,
		NumReadsDiscarded:  discarded,
		NumReadsInProgress: inProgress,
		NumReadsTodo:       len(h.QnameCluster) - (accepted + discarded),
	}

	return g
}
