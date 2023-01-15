MBD2 <- function(D, policy.genton = TRUE )
{
  # D:                 is the dataset containing signals whose depths are to be computed.
  # policy.genton:     specifies whether the computation should follow the genton formula
  #                    or not.

  N = dim(D)[1]
  P = dim(D)[2]

  Depths = rep(0,N);
  if( policy.genton )
  {
    for( i in 1:P )
    {
      rk = rank(D[,i]);

      for( j in 1:N)
      {
        N.up = rk[j]-1;
        N.down = N - rk[j];
        ###################
        # Genton proposal:
        ###################
        Depths[j] = Depths[j] + N.up * N.down + N - 1;
        # Depths are increased by N-1, accounting for the bands given by the current
        # signal and any other signal in the dataset. I prefer the more strict version implemented0
        # hereafter, since in the case of vertically translated signals, the most peripheral
        # would have zero depth, while in the other case their depths would be nonzero (and
        # precisely (N-1)/(N choose 2) )

      }
    }
  }
  else
  {
    for( i in 1:P )
    {
      rk = rank(D[,i]);

      for( j in 1:N)
      {
        N.up = rk[j]-1;
        N.down = N - rk[j];
        Depths[j] = Depths[j] + N.up * N.down;
      }
    }
  }
  Depths = Depths / ( choose(N,2) * P );
  Depths
}

MBD2.rel <- function(D,X, strict.inclusion = FALSE )
{
  # D:                is the reference dataset,
  # X:                is the test set of data whose depths
  #                   have to be computed
  # strict.inclusion: TROVAGLI UN NOME MIGLIORE!!!
  #                   Must the inclusion in the bands be strict or not?

  N = dim(D)[1];
  P = dim(D)[2];

  if( length(dim(as.array(X))) == 1 )
  {
    X = t( as.matrix( X ) )
  }

  N_test = dim(X)[1];

  stopifnot( P == dim(X)[2] );

  Depths = rep(0,N_test);

  if( !strict.inclusion )
  {
    for( i in 1:P )
    {
      for( j in 1 : N_test )
      {
        m.leq = sum( D[,i] <= X[j,i])
        m.geq = sum( D[,i] >= X[j,i])

        exceed = m.leq + m.geq   - N

        if( m.leq < m.geq )
        {
          m.geq = m.geq - exceed
        }
        else
        {
          m.leq = m.leq - exceed
        }

        Depths[j] = Depths[j] + m.leq * m.geq;
      }
    }
  }
  else
  {
    for( i in 1:P )
    {
      for( j in 1 : N_test )
      {
        m.leq = sum( D[,i] < X[j,i])
        m.geq = sum( D[,i] > X[j,i])

        Depths[j] = Depths[j] + m.leq * m.geq;
      }
    }
  }
  Depths = Depths / ( choose(N,2) * P  );
  Depths
}

SD.rel <- function( D, X )
{
  # D:                is the reference dataset,
  # X:                is the test set of data whose depths
  #                   have to be computed


  N = dim(D)[1];
  P = dim(D)[2];

  if( length(dim(as.array(X))) == 1 )
  {
    X = t( as.matrix( X ) )
  }

  N_test = dim(X)[1];

  stopifnot( P == dim(X)[2] );

  Depths = rep(0,N_test);

  for( i.test in 1:N_test )
  {

  	sx = rep(0,P);

#     print(i.test)

  	for( j.ref in 1:N )
  	{
	       nrm = sqrt(sum((X[i.test,] - D[j.ref,])^2));
        if( nrm != 0 )
        {
           sx = sx + (X[i.test,] - D[j.ref,])/nrm;
        }
  	}
  	Depths[i.test] = 1 - norm(sx/N,'2');
  }

  Depths

}
