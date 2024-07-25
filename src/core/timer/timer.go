package timer

import "time"

type Timer struct {
	MapR1                 time.Duration
	MapR2                 time.Duration
	DetermineReadLocation time.Duration
	ExactMatch            time.Duration
}

func NewTimer() *Timer {
	return &Timer{
		MapR1:                 time.Duration(0),
		MapR2:                 time.Duration(0),
		DetermineReadLocation: time.Duration(0),
		ExactMatch:            time.Duration(0),
	}
}
