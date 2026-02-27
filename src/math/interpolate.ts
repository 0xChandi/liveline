import type { LivelinePoint } from '../types'

/**
 * Fritsch-Carlson monotone cubic interpolation at a given time.
 * Uses the same spline math as drawSpline() so the scrub dot
 * sits exactly on the rendered curve.
 * Returns null if time is outside data range.
 */
export function interpolateAtTime(
  points: LivelinePoint[],
  time: number,
): number | null {
  if (points.length === 0) return null
  if (time <= points[0].time) return points[0].value
  if (time >= points[points.length - 1].time) return points[points.length - 1].value

  const n = points.length

  // For only 2 points, linear is correct (drawSpline uses lineTo)
  if (n === 2) {
    const dt = points[1].time - points[0].time
    if (dt === 0) return points[0].value
    const t = (time - points[0].time) / dt
    return points[0].value + (points[1].value - points[0].value) * t
  }

  // 1. Compute secant slopes and intervals (mirrors drawSpline exactly)
  const delta: number[] = new Array(n - 1)
  const h: number[] = new Array(n - 1)
  for (let i = 0; i < n - 1; i++) {
    h[i] = points[i + 1].time - points[i].time
    delta[i] = h[i] === 0 ? 0 : (points[i + 1].value - points[i].value) / h[i]
  }

  // 2. Initial tangent estimates
  const m: number[] = new Array(n)
  m[0] = delta[0]
  m[n - 1] = delta[n - 2]
  for (let i = 1; i < n - 1; i++) {
    if (delta[i - 1] * delta[i] <= 0) {
      m[i] = 0
    } else {
      m[i] = (delta[i - 1] + delta[i]) / 2
    }
  }

  // 3. Fritsch-Carlson constraint: alpha^2 + beta^2 <= 9
  for (let i = 0; i < n - 1; i++) {
    if (delta[i] === 0) {
      m[i] = 0
      m[i + 1] = 0
    } else {
      const alpha = m[i] / delta[i]
      const beta = m[i + 1] / delta[i]
      const s2 = alpha * alpha + beta * beta
      if (s2 > 9) {
        const s = 3 / Math.sqrt(s2)
        m[i] = s * alpha * delta[i]
        m[i + 1] = s * beta * delta[i]
      }
    }
  }

  // 4. Binary search for the interval containing `time`
  let lo = 0
  let hi = n - 1
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1
    if (points[mid].time <= time) lo = mid
    else hi = mid
  }

  // 5. Evaluate cubic Hermite on segment [lo, hi]
  const segH = h[lo]
  if (segH === 0) return points[lo].value
  const t = (time - points[lo].time) / segH
  const t2 = t * t
  const t3 = t2 * t

  // Hermite basis functions
  const h00 = 2 * t3 - 3 * t2 + 1
  const h10 = t3 - 2 * t2 + t
  const h01 = -2 * t3 + 3 * t2
  const h11 = t3 - t2

  return (
    h00 * points[lo].value +
    h10 * segH * m[lo] +
    h01 * points[hi].value +
    h11 * segH * m[hi]
  )
}
