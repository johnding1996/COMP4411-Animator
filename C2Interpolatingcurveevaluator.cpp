#include "C2InterpolatingCurveEvaluator.h"
#include "LinearCurveEvaluator.h"
#include <assert.h>
#include <algorithm>
#include "mat.h"
#include "vec.h"

#define NUM_SEG 20

void C2InterpolatingCurveEvaluator::evaluateCurve(const std::vector<Point>& ptvCtrlPts,
	std::vector<Point>& ptvEvaluatedCurvePts,
	const float& fAniLength,
	const bool& bWrap) const
{
	ptvEvaluatedCurvePts.clear();
	if ((bWrap && ptvCtrlPts.size() < 2) || (!bWrap && ptvCtrlPts.size() < 4))
	{
		LinearCurveEvaluator lineCrvEvl = LinearCurveEvaluator();
		lineCrvEvl.evaluateCurve(ptvCtrlPts, ptvEvaluatedCurvePts, fAniLength, bWrap);
		return;
	}
	const Mat4d ker(
		-1, 3, -3, 1,
		3, -6, 3, 0,
		-3, 3, 0, 0,
		1, 0, 0, 0);
	std::vector<Point> ptvCtrlPtsCpy;
	ptvCtrlPtsCpy.assign(ptvCtrlPts.begin(), ptvCtrlPts.end());
	if (bWrap)
	{
		ptvCtrlPtsCpy.push_back(Point(ptvCtrlPts.front().x + fAniLength, ptvCtrlPts.front().y));
		ptvCtrlPtsCpy.insert(ptvCtrlPtsCpy.begin(), Point(ptvCtrlPts.back().x - fAniLength, ptvCtrlPts.back().y));
	}

	std::vector<Point> ptvCtrlPtsDvr;
	evaluateDerivative(ptvCtrlPtsDvr, ptvCtrlPtsCpy);

	std::vector<Point> ptvBezierCtrlPts;
	for (int i = 0; i < ((int)ptvCtrlPtsCpy.size() - 1); i++)
	{
		ptvBezierCtrlPts.push_back(ptvCtrlPtsCpy[i]);
		ptvBezierCtrlPts.push_back(Point((float)(ptvCtrlPtsCpy[i].x+ ptvCtrlPtsDvr[i].x/3.0),
			(float)(ptvCtrlPtsCpy[i].y + ptvCtrlPtsDvr[i].y / 3.0)));
		ptvBezierCtrlPts.push_back(Point((float)(ptvCtrlPtsCpy[i+1].x - ptvCtrlPtsDvr[i+1].x / 3.0),
			(float)(ptvCtrlPtsCpy[i+1].y - ptvCtrlPtsDvr[i+1].y / 3.0)));
	}
	ptvBezierCtrlPts.push_back(ptvCtrlPtsCpy.back());

	for (int iter = 0; iter <((int)ptvBezierCtrlPts.size() - 3); iter += 3)
	{

		Vec4d ctrl_x(ptvBezierCtrlPts[iter].x, ptvBezierCtrlPts[iter + 1].x,
			ptvBezierCtrlPts[iter + 2].x, ptvBezierCtrlPts[iter + 3].x);
		Vec4d ctrl_y(ptvBezierCtrlPts[iter].y, ptvBezierCtrlPts[iter + 1].y,
			ptvBezierCtrlPts[iter + 2].y, ptvBezierCtrlPts[iter + 3].y);
		ptvEvaluatedCurvePts.push_back(ptvBezierCtrlPts[iter]);
		for (int i = 1; i < NUM_SEG; ++i)
		{
			double t = i / (double)NUM_SEG;
			Vec4d para(t*t*t, t*t, t, 1);
			Point eval_pt((float)(para*ker*ctrl_x), (float)(para*ker*ctrl_y));
			if (eval_pt.x > ptvBezierCtrlPts[iter].x && eval_pt.x < ptvBezierCtrlPts[iter + 3].x
				&& (ptvEvaluatedCurvePts.empty() || eval_pt.x > ptvEvaluatedCurvePts.back().x))
				ptvEvaluatedCurvePts.push_back(eval_pt);
		}
		ptvEvaluatedCurvePts.push_back(ptvBezierCtrlPts[iter + 3]);
	}

	if (!bWrap)
	{
		ptvEvaluatedCurvePts.push_back(Point(0, ptvBezierCtrlPts.front().y));
		ptvEvaluatedCurvePts.push_back(Point(fAniLength, ptvBezierCtrlPts.back().y));
	}
}



void C2InterpolatingCurveEvaluator::evaluateDerivative(std::vector<Point>& ptvCtrlPtsDvr,
	std::vector<Point> ptvCtrlPtsCpy) const
{
	int m = ptvCtrlPtsCpy.size() - 1;
	std::vector<double> coeff(m+1, 0.0);
	std::vector<double> dv_x(m+1, 0.0);
	std::vector<double> dv_y(m+1, 0.0);

	coeff[0] = 0.5;
	for (int i = 1; i < m; ++i) coeff[i] = 1.0 / (4.0 - coeff[i - 1]);
	coeff[m] = 1.0 / (2.0 - coeff[m - 1]);

	dv_y[0] = 1.5 * (ptvCtrlPtsCpy[1].y - ptvCtrlPtsCpy[0].y);
	for (int i = 1; i < m; i++)
		dv_y[i] = coeff[i] * (3 * (ptvCtrlPtsCpy[i + 1].y - ptvCtrlPtsCpy[i - 1].y) - dv_y[i - 1]);
	dv_y[m] = coeff[m] *
		(3 * (ptvCtrlPtsCpy[m].y - ptvCtrlPtsCpy[m - 1].y) - dv_y[m]);

	dv_x[0] = 1.5 * (ptvCtrlPtsCpy[1].x - ptvCtrlPtsCpy[0].x);
	for (int i = 1; i < m; i++)
		dv_x[i] = coeff[i] * (3 * (ptvCtrlPtsCpy[i + 1].x - ptvCtrlPtsCpy[i - 1].x) - dv_x[i - 1]);
	dv_x[m] = coeff[m] *
		(3 * (ptvCtrlPtsCpy[m].x - ptvCtrlPtsCpy[m - 1].x) - dv_x[m]);

	ptvCtrlPtsDvr.push_back(Point((float)(dv_x[m]), (float)(dv_y[m])));
	for (int i = m - 1; i >= 0; i--)
	{
		dv_x[i] = dv_x[i] - coeff[i] * dv_x[i + 1];
		dv_y[i] = dv_y[i] - coeff[i] * dv_y[i + 1];
		ptvCtrlPtsDvr.push_back(Point((float)(dv_x[i]), (float)(dv_y[i])));
	}
	std::reverse(ptvCtrlPtsDvr.begin(), ptvCtrlPtsDvr.end());
}