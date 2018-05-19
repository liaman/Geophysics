//  Copyright 2005-2011 Jason Aaron Osgood, Jean-Francois Poilpret
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package net.java.dev.designgridlayout;

import java.awt.Container;
import java.util.Iterator;
import java.util.List;

import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;

final class SubGrid implements ISubGrid
{
	SubGrid(List<RowItem> items, SubGrid previous, ParentWrapper<Container> parent, 
		JLabel label, int gridspan)
	{
		_items = items;
		_previous = previous;
		_parent = parent;
		_label = label;
		_gridspan = (gridspan <= 0 ? 0 : gridspan);
		if (_label != null)
		{
			_parent.add(_label);
		}
	}

	void spanRow()
	{
		if (_previous == null)
		{
			// Bad use of DesignGridLayout, use a maker component
			add(MarkerHelper.createMarker(1, NO_PREVIOUS_SUBGRID), 1);
		}
		else
		{
			// Find the RowItem, in the above subgrid, matching the current
			// column position.
			RowItem previous = _previous.findItem(_column);
			if (previous == null)
			{
				// Bad use of DesignGridLayout, use a maker component
				add(MarkerHelper.createMarker(1, NO_MATCHING_COMPONENT), 1);
			}
			else
			{
				// It is impossible to say now if this call can succeed, this can
				// only be checked later (checkSpanRows()), for now we consider
				// that it works
				_spanRow = true;
				_items.add(new RowItem(previous));
			}
		}
	}

	void indent(int n)
	{
		_indent = ComponentGapsHelper.instance().getHorizontalIndent() * n;
	}
	
	void add(JComponent child, int span)
	{
		RowItem item;
		if (child != null)
		{
			item = new RowItem(span, child);
			_parent.add(child);
		}
		else
		{
			item = new RowItem(span, EMPTY);
		}
		_column += span;
		_items.add(item);
	}
	
	void checkSpanRows()
	{
		// If there is no remaining spanRow() call to check, then we're done
		if (!_spanRow)
		{
			return;
		}

		// Check that the number of columns in this sub-grid matches the 
		// previous sub-grid. If not, then we have to replace all spanning
		// RowItems with marker components
		if (_previous.gridColumns() != gridColumns())
		{
			for (RowItem item: _items)
			{
				if (!item.isFirstSpanRow())
				{
					JComponent marker = MarkerHelper.createMarker(
						item.span(), UNMATCHED_COLUMNS_SUBGRIDS);
					item.replace(marker);
					_parent.add(marker);
				}
			}
		}
	}
	
	JLabel label()
	{
		return _label;
	}
	
	JComponent leftComponent()
	{
		if (_label != null)
		{
			return _label;
		}
		else
		{
			return (_items.isEmpty() ? null : _items.get(0).component());
		}
	}
	
	@Override public int gridspan()
	{
		return _gridspan;
	}
	
	public void gridspan(int span)
	{
		if (_previous != null && _spanRow)
		{
			_gridspan = _previous.gridspan();
		}
		if (_gridspan == 0)
		{
			_gridspan = span;
		}
	}

	@Override public int labelWidth()
	{
//		return (_label != null ? _label.getPreferredSize().width : 0);
		return (_label != null ? _label.getPreferredSize().width : 0) + _indent;
	}

	@Override public int gridColumns()
	{
		int columns = 0;
		for (RowItem item: _items)
		{
			columns += item.span();
		}
		return columns;
	}

	@Override public int maxColumnWidth(int maxColumns, IExtractor extractor)
	{
		int maxWidth = 0;
		// Note columns (sum item spans), not the count of components
		int columns = gridColumns();
		float divisions = (float) columns / (float) maxColumns;

		for (RowItem item: _items)
		{
			int width = extractor.value(item);

			// Ignores remainder (fudge), which is incorrect if remainder
			// is greater than horizontal gap (hopefully rarely)
			width = (int) ((width * divisions) / item.span());
			maxWidth = Math.max(maxWidth, width);
		}
		return maxWidth;
	}

	public int hgap()
	{
		return ComponentHelper.hgap(_label, _items, _parent.parent());
	}

	public int layoutRow(LayoutHelper helper, int left, int height, int baseline, 
		int hgap, int rowWidth, int labelWidth)
	{
		int x = left;
		int actualHeight = 0;
		// Account for label column
		if (labelWidth > 0)
		{
			if (_label != null)
			{
//				actualHeight = Math.max(0, helper.setSizeLocation(
//					_label, x, labelWidth, height, baseline));
				actualHeight = Math.max(0, helper.setSizeLocation(
					_label, x + _indent, labelWidth - _indent, height, baseline));
			}
			x += labelWidth + hgap;
		}

		int columns = gridColumns();
		if (columns > 0)
		{
			// pre-subtract gaps
			int gridWidth = rowWidth - ((columns - 1) * hgap);
			int columnWidth = gridWidth / columns;

			// fudge is whatever pixels are left over
			int fudge = gridWidth % columns;

			Iterator<RowItem> i = _items.iterator();
			while (i.hasNext())
			{
				RowItem item = i.next();
				int span = item.span();
				int width = columnWidth * span + ((span - 1) * hgap);
				// Apply the fudge to the last component/column
				if (!i.hasNext())
				{
					width += fudge;
				}
				if (item.isFirstSpanRow())
				{
					JComponent component = item.component();
					actualHeight = Math.max(0, helper.setSizeLocation(
						component, x, width, height, baseline));
				}
				x += width + hgap;
			}
		}
		return actualHeight;
	}

	public List<RowItem> items()
	{
		return _items;
	}
	
	private RowItem findItem(int column)
	{
		int i = 0;
		for (RowItem item: _items)
		{
			if (i == column)
			{
				return item;
			}
			i += item.span();
			if (i > column)
			{
				return null;
			}
		}
		return null;
	}
	
	static final private String NO_PREVIOUS_SUBGRID = 
		"spanRow() cannot work on a grid-row with no grid-row immediately above, " +
		"or with no matching sub-grid (same column position) in the above grid-row";
	static final private String NO_MATCHING_COMPONENT = 
		"spanRow() cannot work when there is no component, on the above grid-row, " +
		"with a matching column location";
	static final private String UNMATCHED_COLUMNS_SUBGRIDS = 
		"spanRow() cannot work on a sub-grid where the number of columns is different " +
		"from the above sub-grid";

	static final private JComponent EMPTY = new JPanel();
	static
	{
		EMPTY.setOpaque(false);
	}

	final private List<RowItem> _items;
	final private SubGrid _previous;
	final private ParentWrapper<Container> _parent;
	final private JLabel _label;
	// 0 means auto span until the right-most edge
	private int _gridspan;
	private boolean _spanRow;
	private int _column = 0;
	private int _indent = 0;
}
