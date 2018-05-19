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

// Fake subgrid that is used whenever a real SubGrid spans several grids
final class EmptySubGrid implements ISubGrid
{
	@Override public int gridColumns()
	{
		return 0;
	}

	@Override public int labelWidth()
	{
		return 0;
	}

	@Override public int maxColumnWidth(int maxColumns, IExtractor extractor)
	{
		return 0;
	}

	@Override public int gridspan()
	{
		return 1;
	}
}
