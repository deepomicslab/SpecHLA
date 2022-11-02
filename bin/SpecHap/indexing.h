//
// Created by yyh on 11/27/2018.
// credit to document for Eigen 3.3.5: https://eigen.tuxfamily.org/dox/TopicCustomizing_NullaryExpr.html
//

#ifndef PHASER_INDEXING_H
#define PHASER_INDEXING_H
#include "Eigen/Dense"

template<class ArgType, class RowIndexType, class ColIndexType>
class indexing_functor {
    const ArgType &m_arg;
    const RowIndexType &m_rowIndices;
    const ColIndexType &m_colIndices;
public:
    typedef Eigen::Matrix<typename ArgType::Scalar,
            RowIndexType::SizeAtCompileTime,
            ColIndexType::SizeAtCompileTime,
            ArgType::Flags&Eigen::RowMajorBit?Eigen::RowMajor:Eigen::ColMajor,
            RowIndexType::MaxSizeAtCompileTime,
            ColIndexType::MaxSizeAtCompileTime> MatrixType;
    indexing_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
            : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
    {}
    const typename ArgType::Scalar& operator() (Eigen::Index row, Eigen::Index col) const
    {
        return m_arg(m_rowIndices[row], m_colIndices[col]);
    }
};


template <class ArgType, class RowIndexType, class ColIndexType>
Eigen::CwiseNullaryOp<indexing_functor<ArgType,RowIndexType,ColIndexType>, typename indexing_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
mat_indexing(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
    typedef indexing_functor<ArgType, RowIndexType, ColIndexType> Func;
    typedef typename Func::MatrixType MatrixType;
    return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(),
                                   Func(arg.derived(), row_indices, col_indices));
}

#endif //PHASER_INDEXING_H
